nextflow.enable.dsl = 2

include { MUTATION_CALLING }  from './subworkflows/mutation_calling'
include { CONSENSUS_CALLING } from './subworkflows/consensus_calling'

workflow {

  def manifest = params.pairs_tsv ?: params.samplesheet
  if( !manifest ) error "Missing --pairs_tsv (legacy alias: --samplesheet)"

  def inferBaiPath = { String bamPath ->
    def candidates = []
    if( bamPath?.endsWith('.bam') ) {
      candidates << bamPath.replaceFirst(/\.bam$/, '.bai')
    }
    candidates << "${bamPath}.bai"
    def existing = candidates.find { new File(it).exists() }
    return existing ?: candidates[0]
  }

  ch_samples = Channel
    .fromPath(manifest)
    .splitCsv(header:true, sep: manifest.toString().toLowerCase().endsWith('.csv') ? ',' : '\t')
    .map { row ->
      // Normalize header/value whitespace and strip optional UTF-8 BOM on first header.
      def norm = row.collectEntries { k, v ->
        def key = k?.toString()?.replaceFirst('^\\uFEFF', '')?.trim()
        def val = v?.toString()?.trim()
        [(key): val]
      }

      if( !norm.pair_id || !norm.tumor_id || !norm.tumor_bam || !norm.normal_id || !norm.normal_bam || !norm.gender ) {
        error "Pair manifest must include: pair_id,tumor_id,tumor_bam,normal_id,normal_bam,gender"
      }

      def tumorBam = norm.tumor_bam.toString()
      def normalBam = norm.normal_bam.toString()
      def tumorBai = inferBaiPath(tumorBam)
      def normalBai = inferBaiPath(normalBam)

      tuple(
        norm.tumor_id.toString(),
        file(tumorBam),
        file(tumorBai),
        file(normalBam),
        file(normalBai),
        norm.normal_id.toString()
      )
    }

  ch_called = MUTATION_CALLING(ch_samples)
  CONSENSUS_CALLING(ch_called)
}
