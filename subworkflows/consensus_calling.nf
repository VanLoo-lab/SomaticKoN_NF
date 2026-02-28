nextflow.enable.dsl = 2

include {
  NORMALIZE_VCF as NORMALIZE_VCF_M2
  NORMALIZE_VCF as NORMALIZE_VCF_MS
  NORMALIZE_VCF as NORMALIZE_VCF_ST
  SPLIT_MNV as SPLIT_MNV_M2
  SPLIT_MNV as SPLIT_MNV_MS
  SPLIT_MNV as SPLIT_MNV_ST
  STANDARDIZE_VCF as STANDARDIZE_VCF_M2
  STANDARDIZE_VCF as STANDARDIZE_VCF_MS
  STANDARDIZE_VCF as STANDARDIZE_VCF_ST
  LABEL_EVENTS as LABEL_EVENTS_M2
  LABEL_EVENTS as LABEL_EVENTS_MS
  LABEL_EVENTS as LABEL_EVENTS_ST
  CONSENSUS_K_OF_N
  MERGE_EVENTS
  MAKE_SOMATIC_METRICS
  PLOT_SOMATIC_QC
  RENDER_QUARTO
} from '../workflows/consensus_processes'

workflow CONSENSUS_CALLING {

  take:
    ch_called   // tuple(sample, normal_id, mutect2_raw_vcf, muse2_raw_vcf, strelka2_raw_vcf)

  main:

    // ---- per-caller postprocess: normalize -> split -> std -> label
    ch_m2_evt = ch_called
      .map { s, normal_id, m2, ms, st -> tuple(s, normal_id, 'mutect2', m2) }
      | NORMALIZE_VCF_M2
      | SPLIT_MNV_M2
      | STANDARDIZE_VCF_M2
      | LABEL_EVENTS_M2

    ch_ms_evt = ch_called
      .map { s, normal_id, m2, ms, st -> tuple(s, normal_id, 'muse2', ms) }
      | NORMALIZE_VCF_MS
      | SPLIT_MNV_MS
      | STANDARDIZE_VCF_MS
      | LABEL_EVENTS_MS

    ch_st_evt = ch_called
      .map { s, normal_id, m2, ms, st -> tuple(s, normal_id, 'strelka2', st) }
      | NORMALIZE_VCF_ST
      | SPLIT_MNV_ST
      | STANDARDIZE_VCF_ST
      | LABEL_EVENTS_ST

    // ---- consensus (needs 3 labeled event vcfs per sample)
    ch_cons = ch_m2_evt
      .join(ch_ms_evt, by: 0)
      .join(ch_st_evt, by: 0)
      .map { s, m2evt, msevt, stevt -> tuple(s, m2evt, msevt, stevt) }
      | CONSENSUS_K_OF_N

    ch_merged = ch_cons | MERGE_EVENTS

    // ---- metrics needs raw vcfs + consensus + merged
    ch_metrics_input = ch_called
      .join(ch_cons.map { s, cons_vcf, kpass -> tuple(s, cons_vcf) }, by: 0)
      .join(ch_merged.map { s, merged_vcf, kpass -> tuple(s, merged_vcf, kpass) }, by: 0)
      .map { s, normal_id, m2raw, msraw, straw, cons_vcf, merged_vcf, kpass ->
        tuple(s, m2raw, msraw, straw, cons_vcf, merged_vcf, kpass)
      }

    ch_metrics = MAKE_SOMATIC_METRICS(ch_metrics_input)
    ch_plots   = PLOT_SOMATIC_QC(ch_metrics)
    ch_html    = RENDER_QUARTO(ch_plots)

  emit:
    ch_cons
    ch_merged
    ch_html
}
