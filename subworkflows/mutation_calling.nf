nextflow.enable.dsl = 2

include {
  CHECK_BAI
  MUTECT2_SCATTER
  MUTECT2_GATHER
  MUSE2_CALL
  STRELKA2_CALL
} from '../workflows/mutation_processes'

workflow MUTATION_CALLING {

  take:
    ch_samples   // tuple(sample, tumor_bam, tumor_bai, normal_bam, normal_bai, normal_id?)

  main:
    // Create chromosome channel for scattering
    ch_chroms = Channel.fromList(params.chroms)

    // 0) Validate BAM index presence (fails fast with a helpful error)
    ch_checked = CHECK_BAI(ch_samples)

    // 1) Mutect2 scatter -> gather
    ch_scatter = MUTECT2_SCATTER(ch_checked, ch_chroms)

    ch_scatter_grouped = ch_scatter
      .groupTuple(by: 0)
      .map { sample, chrs, scatter_files ->
        tuple(sample, scatter_files.flatten())
      }

    ch_mutect2 = ch_checked
      .join(ch_scatter_grouped, by: 0)
      .map { sample, tumor_bam, tumor_bai, normal_bam, normal_bai, normal_id, scatter_files ->
        tuple(sample, tumor_bam, tumor_bai, normal_bam, normal_bai, normal_id, scatter_files)
      }
      | MUTECT2_GATHER

    // 2) MuSE2 + Strelka2 run in parallel
    ch_muse2    = MUSE2_CALL(ch_checked)
    ch_strelka2 = STRELKA2_CALL(ch_checked)

    // 3) Join raw VCFs by sample
    ch_join = ch_mutect2
      .join(ch_muse2, by: 0)
      .join(ch_strelka2, by: 0)
      .map { s, mutect2_vcf, muse2_vcf, strelka2_vcf ->
        tuple(s, mutect2_vcf, muse2_vcf, strelka2_vcf)
      }

    // 4) Attach normal_id for downstream standardization
    ch_join_with_normal = ch_join
      .join(ch_checked.map { sample, tumor_bam, tumor_bai, normal_bam, normal_bai, normal_id -> tuple(sample, normal_id) }, by: 0)
      .map { s, mutect2_vcf, muse2_vcf, strelka2_vcf, normal_id ->
        tuple(s, normal_id, mutect2_vcf, muse2_vcf, strelka2_vcf)
      }

  emit:
    ch_join_with_normal
}
