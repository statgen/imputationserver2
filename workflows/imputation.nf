nextflow.enable.types = true

include { RefPanelParams ; ImputationParams } from '../types/params'
include { VcfMetafile ; MinimacInput ; ImputedChunk } from '../types/workflow'

include { MINIMAC4 } from '../modules/local/imputation/minimac4'

workflow IMPUTATION {
    take:
    refpanel: RefPanelParams
    imputation: ImputationParams
    phased_ch: Channel<VcfMetafile>
    threads: Integer

    main:
    def minimac_map: Path? = (refpanel.mapMinimac == null)
        ? null
        : file(refpanel.mapMinimac, checkIfExists: true)

    def chromosomes = channel.fromList(1..22)
        .map { chr -> "${chr}" }
        .mix(channel.of('X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT'))

    def minimac_m3vcf_ch = chromosomes
        .map { chr ->
            record(
                chr: chr,
                m3vcf: file(PatternUtil.parse(refpanel.genotypes, [chr: chr])),
            )
        }
        .filter { it -> it.m3vcf.exists() }

    def phased_m3vcf_ch: Channel<MinimacInput> = phased_ch.join(minimac_m3vcf_ch, by: 'chr')

    minimac = MINIMAC4(
        phased_m3vcf_ch,
        minimac_map,
        refpanel.build,
        imputation.window,
        imputation.minimac_min_ratio,
        imputation.min_r2,
        imputation.decay,
        imputation.diff_threshold,
        imputation.prob_threshold,
        imputation.prob_threshold_s1,
        imputation.min_recom,
        threads,
    )

    emit:
    imputed_chunks: Channel<ImputedChunk> = minimac
}
