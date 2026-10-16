nextflow.enable.types = true

include { RefPanelParams ; ImputationParams ; PhasingParams } from '../types/params'
include { INPUT_VALIDATION_VCF } from '../modules/local/input_validation/input_validation_vcf'

workflow INPUT_VALIDATION {
    take:
    population: String
    files: Channel<Path>
    refpanel: RefPanelParams
    imputation: ImputationParams
    phasing: PhasingParams
    build: String
    mode: String
    chunksize: Integer
    min_samples: Integer
    max_samples: Integer

    main:
    result = INPUT_VALIDATION_VCF(
        population,
        files.collect(),
        refpanel,
        imputation,
        phasing,
        build,
        mode,
        chunksize,
        min_samples,
        max_samples,
    )

    emit:
    validated_files = result.validated_files
    validation_report = result.validation_report
}
