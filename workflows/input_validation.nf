include { INPUT_VALIDATION_VCF } from '../modules/local/input_validation/input_validation_vcf'

workflow INPUT_VALIDATION {
    take:
    population // String
    files // Channel<Path>
    refpanel//: RefPanelParams

    main:
    INPUT_VALIDATION_VCF(population, files.collect(), refpanel)

    emit:
    validated_files = INPUT_VALIDATION_VCF.out.validated_files.collect()
    validation_report = INPUT_VALIDATION_VCF.out.validation_report
}
