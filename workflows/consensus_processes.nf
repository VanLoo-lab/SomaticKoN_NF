nextflow.enable.dsl = 2

def loadRefs = """
source "${params.resolveRefs}" --ref "${params.ref}" --refdir "${params.refdir}" --seq "${params.seq}"
"""

process NORMALIZE_VCF {
  tag { sample }

  cpus { params.threads }
  memory '100 GB'
  time '6h'

  input:
    tuple val(sample), val(normal_id), val(caller), path(raw_vcf)

  output:
    tuple val(sample), val(normal_id), val(caller), path("${sample}.${caller}.norm.vcf.gz"), path("${sample}.${caller}.presplit.vcf.gz")

  script:
    """
    source "${params.libDir}/conda_utils.sh"
  
    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    ${loadRefs}
    bash "${params.binDir}/normalize_vcf.sh" \
      --in-vcf "${raw_vcf}" \
      --ref "\$REF_GENOME" \
      --out "${sample}.${caller}.norm.vcf.gz" \
      --out-presplit "${sample}.${caller}.presplit.vcf.gz" \
      --threads ${params.threads}
    """
}

process SPLIT_MNV {
  tag { sample }

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), val(normal_id), val(caller), path(norm_vcf), path(presplit_vcf)

  output:
    tuple val(sample), val(normal_id), val(caller), path("${sample}.${caller}.split.vcf.gz"), path(presplit_vcf)

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    bash "${params.binDir}/split_mnv.sh" \
      --in-vcf "${norm_vcf}" \
      --out "${sample}.${caller}.split.vcf.gz" \
      --workdir .
    """
}

process STANDARDIZE_VCF {
  tag { sample }

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), val(normal_id), val(caller), path(split_vcf), path(presplit_vcf)

  output:
    tuple val(sample), val(normal_id), val(caller), path("${sample}.${caller}.std.vcf.gz"), path(presplit_vcf)

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    bash "${params.binDir}/standardize_vcf_samples.sh" \
      --in-vcf "${split_vcf}" \
      --out-vcf "${sample}.${caller}.std.vcf.gz" \
      --tumor-hint "${sample}" \
      --normal-hint "${normal_id ?: (params.normal_id ?: 'NORMAL')}" \
      --tumor-name TUMOR \
      --normal-name NORMAL
    """
}

process LABEL_EVENTS {
  tag { sample }

  cpus { params.threads }
  memory '100 GB'
  time '4h'

  input:
    tuple val(sample), val(normal_id), val(caller), path(std_vcf), path(presplit_vcf)

  output:
    tuple val(sample), path("${sample}.${caller}.event.vcf.gz")

  script:
    """
    source "${params.libDir}/conda_utils.sh"
    
    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    bash "${params.binDir}/label_events_from_presplit.sh" \
      --presplit "${presplit_vcf}" \
      --final "${std_vcf}" \
      --out "${sample}.${caller}.event.vcf.gz" \
      --threads ${params.threads}
    """
}

process CONSENSUS_K_OF_N {
  tag { sample }
  publishDir "${params.outdir}/${sample}/consensus", mode: 'copy', overwrite: true

  cpus { params.threads }
  memory '100 GB'
  time '8h'

  input:
    tuple val(sample), path(m2_evt), path(ms_evt), path(st_evt)

  output:
    tuple val(sample), path("${sample}.kn.k${params.k_pick}.vcf"), path("kpass")

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    ${loadRefs}

    workdir="${sample}.kn.work"
    mkdir -p "\${workdir}" kpass

    bash "${params.binDir}/consensus_k_of_n.sh" \
      --mode compare \
      --ref "\$REF_GENOME" \
      --workdir "\${workdir}" \
      --vcf mutect2="${m2_evt}" \
      --vcf muse2="${ms_evt}" \
      --vcf strelka2="${st_evt}" \
      --out-prefix "\${workdir}/${sample}.kn" \
      --k-range "${params.k_range}"

    test -s "\${workdir}/${sample}.kn.k${params.k_pick}.vcf"
    cp "\${workdir}/${sample}.kn.k${params.k_pick}.vcf" "${sample}.kn.k${params.k_pick}.vcf"

    shopt -s nullglob
    for f in "\${workdir}"/k*.pass.tsv; do cp "\$f" "kpass/"; done
    """
}

process MERGE_EVENTS {
  tag { sample }
  publishDir "${params.outdir}/${sample}/consensus", mode: 'copy', overwrite: true

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), path(cons_vcf), path(kpass_dir)

  output:
    tuple val(sample), path("${sample}.kn.merged.vcf"), path(kpass_dir)

  script:
    """
    source "${params.libDir}/conda_utils.sh"
    
    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    workdir="${sample}.merge.work"
    mkdir -p "\${workdir}"

    bash "${params.binDir}/merge_events.sh" \
      --in "${cons_vcf}" \
      --out "${sample}.kn.merged.vcf" \
      --workdir "\${workdir}"

    test -s "${sample}.kn.merged.vcf"
    """
}

process MAKE_SOMATIC_METRICS {
  tag { sample }

  // Do not publish intermediate reports/metrics; only final report folder.

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), path(m2_raw), path(ms_raw), path(st_raw), path(cons_vcf), path(merged_vcf), path(kpass_dir)

  output:
    tuple val(sample), path("reports")

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    mkdir -p reports/metrics

    shopt -s nullglob
    kpass_args=()
    for f in "${kpass_dir}"/k*.pass.tsv; do
      fname="\$(basename "\$f")"
      k="\${fname#k}"
      k="\${k%.pass.tsv}"
      if [[ "\${k}" =~ ^[0-9]+\$ ]]; then
        kpass_args+=(--kpass "k=\${k} path=\$f")
      fi
    done

    bash "${params.binDir}/make_somatic_metrics.sh" \
      --outdir reports/metrics \
      --sample-id "${sample}" \
      --vcf "stage=raw caller=mutect2 path=${m2_raw}" \
      --vcf "stage=raw caller=muse2   path=${ms_raw}" \
      --vcf "stage=raw caller=strelka2 path=${st_raw}" \
      --vcf "stage=consensus caller=2of3 path=${cons_vcf}" \
      --vcf "stage=merged caller=2of3 path=${merged_vcf}" \
      --anchor-vcf "stage=raw caller=mutect2 path=${m2_raw}" \
      --dp-sample "${sample}" \
      "\${kpass_args[@]}"
    """
}

process PLOT_SOMATIC_QC {
  tag { sample }
  
  // Do not publish intermediate plots; only final report folder.

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), path(report_dir)

  output:
    tuple val(sample), path(report_dir)

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs-plotting \
    "${params.envDir}/plotting.yml" \
    "${params.condaEnvDir}"

    mkdir -p "${report_dir}/plots"

    Rscript "${params.binDir}/plot_somatic_qc.R" \
      "${report_dir}/metrics" \
      "${report_dir}/plots" \
      --sample-id "${sample}" \
      --format png
    """
}

process RENDER_QUARTO {
  tag { sample }
  publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true, saveAs: { name ->
    if (name == 'site') return 'final_report'
    return null
  }

  cpus 1
  memory '100 GB'
  time '2h'

  input:
    tuple val(sample), path(report_dir)

  output:
    path("site")

  script:
    """
    source "${params.libDir}/conda_utils.sh"

    ensure_conda_env \
    wgs-plotting \
    "${params.envDir}/plotting.yml" \
    "${params.condaEnvDir}"

    workdir="\$PWD"
    cd "${report_dir}"

    mkdir -p "\${workdir}/site"

    cp "${params.binDir}/somatic_report.qmd" ./somatic_report.qmd
    cp "${params.binDir}/custom.css" ./custom.css

    quarto render somatic_report.qmd \
      -P sample_id:"${sample}" \
      --output-dir "\${workdir}/site" \
      -o somatic_report.html

    # site is already in the task workdir for publishing
    """
}
