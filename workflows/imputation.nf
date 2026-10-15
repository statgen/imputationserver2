// nextflow.enable.types =  true

// include { RefPanelParams; ImputationParams } from '../types/params'
include { MINIMAC4 } from '../modules/local/imputation/minimac4'

workflow IMPUTATION {
    take:
    refpanel//: RefPanelParams
    imputation//: ImputationParams
    phased_ch

    main:
    if (refpanel.mapMinimac == null) {
        minimac_map = []
    } else {
        minimac_map = file(refpanel.mapMinimac, checkIfExists: true)
    }

    chromosomes = channel.of(1..22, 'X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT')
    minimac_m3vcf_ch = chromosomes
        .map {
            it ->
                def genotypes_file = file(PatternUtil.parse(refpanel.genotypes, [chr: it]))
                    if(!genotypes_file.exists()){
                        return null
                    }
                return tuple(it.toString(), genotypes_file)
        }

    phased_m3vcf_ch = phased_ch.combine(minimac_m3vcf_ch, by: 0)

    MINIMAC4(
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
    )

    emit:
    imputed_chunks = MINIMAC4.out.imputed_chunks
}
