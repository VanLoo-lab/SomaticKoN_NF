# WGS Somatic Variant Calling Pipeline (Nextflow DSL2)

This repository contains a Nextflow DSL2 pipeline for tumor-normal somatic variant calling and consensus reporting.

## Overview

The pipeline runs three callers in parallel:

- Mutect2 (scatter by chromosome, then gather)
- MuSE2
- Strelka2

Then it performs a consensus workflow:

- Normalize/split/standardize/label per caller
- Build K-of-N consensus
- Merge event-level variants
- Generate QC metrics and plots
- Render a final Quarto HTML report

Main entrypoint:

- `main.nf`

Configuration:

- `conf/nextflow.config`

## Workflow

Top-level flow:

1. Parse pair manifest (`pair_id,tumor_id,tumor_bam,normal_id,normal_bam,gender`)
2. `MUTATION_CALLING` subworkflow
3. `CONSENSUS_CALLING` subworkflow

Subworkflows:

- `subworkflows/mutation_calling.nf`
- `subworkflows/consensus_calling.nf`

Process definitions:

- `workflows/mutation_processes.nf`
- `workflows/consensus_processes.nf`

## Requirements

- Nextflow with DSL2 support
- LSF cluster (default executor in config)
- Conda or Mamba (conda is enabled in config)
- Reference bundle expected by `bin/lib/resolve_refs.sh`

Environment files:

- `envs/mutect2.yml`
- `envs/muse2.yml`
- `envs/strelka2.yml`
- `envs/plotting.yml`

## Installation

1. Clone the repository and enter it.

```bash
git clone <repo-url>
cd pipelines
```

2. Install Nextflow (choose one method appropriate for your environment).

```bash
# Example: with Conda
conda create -n nextflow -c conda-forge -c bioconda nextflow
conda activate nextflow
```

3. Create the pipeline Conda environments in a dedicated env root directory.

```bash
# Choose an env root directory (example)
mkdir -p /path/to/conda_envs

conda env create --prefix /path/to/conda_envs/wgs          --file envs/mutect2.yml
conda env create --prefix /path/to/conda_envs/muse2        --file envs/muse2.yml
conda env create --prefix /path/to/conda_envs/strelka2     --file envs/strelka2.yml
conda env create --prefix /path/to/conda_envs/wgs-plotting --file envs/plotting.yml
```

4. Add the Conda environment directory to Nextflow so the pipeline can find these envs.

Option A (recommended): set it in `conf/nextflow.config`.

```groovy
params {
  condaEnvDir = "/path/to/conda_envs"
}
```

Option B: pass it at runtime.

```bash
nextflow run main.nf -c conf/nextflow.config --condaEnvDir /path/to/conda_envs ...
```

Notes:

- `params.condaEnvDir` is optional, but setting it explicitly avoids ambiguity.
- If not set, the pipeline probes: `${params.pipeline}/conda_envs`, `$HOME/.conda/envs`, and `$HOME/conda-envs`.

## Input

### Pair Manifest

Provide `--pairs_tsv` as TSV (or CSV). Required columns:

- `pair_id`
- `tumor_id`
- `tumor_bam`
- `normal_id`
- `normal_bam`
- `gender`

Notes:

- `tumor_bai` / `normal_bai` columns are no longer required.
- BAI paths are inferred from BAM paths (`<bam>.bai` or `<bam-without-.bam>.bai`).
- Legacy `--samplesheet` is still accepted as an alias to `--pairs_tsv`.

Example TSV:

```tsv
pair_id	tumor_id	tumor_bam	normal_id	normal_bam	gender
MK545-A_pair	MK545-A	/path/MK545-A.bam	MK545-Control	/path/MK545-Control.bam	male
```

## Running

### Local profile

```bash
nextflow run main.nf \
  -profile local \
  -c conf/nextflow.config \
  --pairs_tsv tests/samplesheet.tsv \
  --outdir tests/nf_out
```

### LSF profile

```bash
nextflow run main.nf \
  -profile lsf \
  -c conf/nextflow.config \
  --pairs_tsv tests/samplesheet.tsv \
  --outdir tests/nf_out
```

### Resume

```bash
nextflow run main.nf -c conf/nextflow.config -resume ...
```

Notes:

- Resume requires the same launch/work context and unchanged task signatures.
- Existing trace/report/timeline/dag files are configured to overwrite in `conf/nextflow.config`.

## Key Parameters

From `conf/nextflow.config`:

- `params.pairs_tsv`: pair manifest path
- `params.samplesheet`: legacy alias for `pairs_tsv`
- `params.outdir`: output directory root
- `params.ref`: reference key (default `hg38`)
- `params.refdir`: reference base directory
- `params.seq`: sequencing type (default `WGS`)
- `params.pipeline`: pipeline root directory (default `projectDir`)
- `params.binDir`: tools scripts path (`${params.pipeline}/modules/somatic`)
- `params.resolveRefs`: reference resolver script (`${params.pipeline}/bin/lib/resolve_refs.sh`)
- `params.libDir`: shared library scripts (`${params.pipeline}/bin/lib`)
- `params.envDir`: env YAML directory (`${params.pipeline}/envs`)
- `params.threads`: default thread count
- `params.chroms`: scatter chromosomes (default `chr1-22, chrX`)
- `params.k_range`: consensus K range
- `params.k_pick`: selected K for final consensus
- `params.condaEnvDir`: conda envs location

## Output Layout

Published outputs are organized by sample under `--outdir`:

- `<outdir>/<sample>/mutect2/scatter/`
- `<outdir>/<sample>/mutect2/`
- `<outdir>/<sample>/muse2/`
- `<outdir>/<sample>/strelka2/`
- `<outdir>/<sample>/consensus/`
- `<outdir>/<sample>/final_report/`

The final HTML report is rendered and published under:

- `<outdir>/<sample>/final_report/somatic_report.html`

Pipeline run metadata is written to:

- `<outdir>/pipeline_info/timeline.html`
- `<outdir>/pipeline_info/report.html`
- `<outdir>/pipeline_info/trace.tsv`
- `<outdir>/pipeline_info/dag.svg`

## Execution Behavior

- Default executor: `lsf`
- Queue: `cgel`
- `executor.perJobMemLimit = true` is enabled (LSF memory handled per job)
- Caller processes (`MUTECT2_SCATTER`, `MUTECT2_GATHER`, `MUSE2_CALL`, `STRELKA2_CALL`) use `errorStrategy = 'ignore'` so one caller failure does not immediately terminate the other caller branches

## Troubleshooting

- Missing output errors usually indicate a mismatch between `output:` pattern and actual script filenames.
- If resume appears to rerun everything, ensure:
  - Same working directory context
  - Same config/parameters where possible
  - Existing `work/` and `.nextflow/` cache accessible
- Do not delete `work/` if you need resume/debug support.

## Repository Structure

- `main.nf`: entry workflow
- `conf/nextflow.config`: parameters, executors, defaults
- `workflows/`: process definitions
- `subworkflows/`: orchestration of process groups
- `modules/somatic/`: caller and utility shell scripts
- `bin/lib/`: shared helper scripts (including reference resolver and conda helpers)
- `envs/`: conda environment specs
- `tests/`: test data, run scripts, and local run artifacts
