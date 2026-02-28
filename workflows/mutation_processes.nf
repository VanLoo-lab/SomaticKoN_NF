nextflow.enable.dsl = 2

def loadRefs = """
source "${params.resolveRefs}" --ref "${params.ref}" --refdir "${params.refdir}" --seq "${params.seq}"
"""

process CHECK_BAI {
  tag { sample }

  // no conda needed
  cpus 1
  memory '100 GB'
  time '10m'

  input:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id)

  output:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id)

  script:
    """
    set -euo pipefail

    if [[ ! -s "${tumor_bai}" ]]; then
      echo "ERROR: Missing BAI for tumor BAM: ${tumor_bam}" >&2
      echo "Provided tumor BAI: ${tumor_bai}" >&2
      exit 2
    fi

    if [[ ! -s "${normal_bai}" ]]; then
      echo "ERROR: Missing BAI for normal BAM: ${normal_bam}" >&2
      echo "Provided normal BAI: ${normal_bai}" >&2
      exit 2
    fi
    """
}

process MUTECT2_SCATTER {
  tag { "${sample}:${chr}" }
  publishDir "${params.outdir}/${sample}/mutect2/scatter", mode: 'copy', overwrite: true
  
  cpus 1
  memory '100 GB'
  time '20h'

  input:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id)
    each chr

  output:
    tuple val(sample), val(chr), path("${sample}.mutect2.${chr}.*")

  script:
    """
    source "${params.libDir}/conda_utils.sh"
  
    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    ${loadRefs}

    bash "${params.binDir}/mutect2_scatter.sh" \
      --sample-id "${sample}" \
      --normal-id "${normal_id ?: 'NORMAL'}" \
      --tumor-bam "${tumor_bam}" \
      --normal-bam "${normal_bam}" \
      --ref "\$REF_GENOME" \
      --interval "${chr}" \
      --germline-resource "\$GNOMAD" \
      --pon "\$PON" \
      --out-prefix "${sample}.mutect2.${chr}"
    """
}

process MUTECT2_GATHER {
  tag { sample }
  publishDir "${params.outdir}/${sample}/mutect2", mode: 'copy', overwrite: true

  cpus 1
  memory '200 GB'
  time '10h'

  input:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id), path(scatter_files)

  output:
    tuple val(sample), path("${sample}.mutect2.pass.vcf.gz")

  script:
    """
    ${loadRefs}

    source "${params.libDir}/conda_utils.sh"
  
    ensure_conda_env \
    wgs \
    "${params.envDir}/mutect2.yml" \
    "${params.condaEnvDir}"

    all_f1r2_input=\$(for c in ${params.chroms.join(' ')}; do printf -- "-I ${sample}.mutect2.\${c}.f1r2.tar.gz "; done)
    all_vcf_input=\$(for c in ${params.chroms.join(' ')}; do printf -- "${sample}.mutect2.\${c}.unfiltered.vcf.gz "; done)
    all_stat_input=\$(for c in ${params.chroms.join(' ')}; do printf -- "-stats ${sample}.mutect2.\${c}.unfiltered.vcf.gz.stats "; done)

    bash "${params.binDir}/mutect2_gather.sh" \
      --sample-id "${sample}" \
      --normal-id "${normal_id ?: 'NORMAL'}" \
      --ref "\$REF_GENOME" \
      --tumor-bam "${tumor_bam}" \
      --normal-bam "${normal_bam}" \
      --common-snps "\$COMMONSNPS" \
      --scatter-vcfs "\${all_vcf_input}" \
      --scatter-stats "\${all_stat_input}" \
      --scatter-f1r2 "\${all_f1r2_input}" \
      --out-prefix "${sample}.mutect2" \
      --pass-only true

    test -s "${sample}.mutect2.pass.vcf.gz"
    """
}

process MUSE2_CALL {
  tag { sample }
  publishDir "${params.outdir}/${sample}/muse2", mode: 'copy', overwrite: true

  cpus { params.threads }
  memory '200 GB'
  time '20h'

  input:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id)

  output:
    tuple val(sample), path("${sample}.muse2.pass.vcf.gz")

  script:
    """
    source "${params.libDir}/conda_utils.sh"
  
    ensure_conda_env \
    muse2 \
    "${params.envDir}/muse2.yml" \
    "${params.condaEnvDir}"

    ${loadRefs}

    bash "${params.binDir}/muse2_call.sh" \
      --sample-id "${sample}" \
      --tumor-bam "${tumor_bam}" \
      --normal-bam "${normal_bam}" \
      --ref "\$REF_GENOME" \
      --dbsnp "\$DBSNP" \
      --outdir . \
      --threads ${params.threads} \
      --exome false \
      --keep-all false

    test -s "${sample}.muse2.pass.vcf.gz"
    """
}


process STRELKA2_CALL {
  tag { sample }
  publishDir "${params.outdir}/${sample}/strelka2", mode: 'copy', overwrite: true

  cpus { params.threads }
  memory '200 GB'
  time '20h'

  input:
    tuple val(sample), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(normal_id)

  output:
    tuple val(sample), path("${sample}.strelka2.pass.vcf.gz")

  script:
    """
    source "${params.libDir}/conda_utils.sh"
  
    ensure_conda_env \
    strelka2 \
    "${params.envDir}/strelka2.yml" \
    "${params.condaEnvDir}"

    ${loadRefs}
    mkdir -p strelka2

    bash "${params.binDir}/strelka2_call.sh" \
      --sample-id "${sample}" \
      --tumor-bam "${tumor_bam}" \
      --normal-bam "${normal_bam}" \
      --ref "\$REF_GENOME" \
      --outdir "strelka2" \
      --threads ${params.threads} \
      --exome false \
      --keep-all false

    if [ -s "strelka2/${sample}.strelka2.pass.vcf.gz" ]; then
      cp "strelka2/${sample}.strelka2.pass.vcf.gz" "${sample}.strelka2.pass.vcf.gz"
    else
      cp "strelka2/${sample}/${sample}.strelka2.pass.vcf.gz" "${sample}.strelka2.pass.vcf.gz"
    fi
    """
}
