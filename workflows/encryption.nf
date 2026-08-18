include { COMPRESSION_ENCRYPTION_VCF } from '../modules/local/compression/compression_encryption_vcf'
include { COMPRESSION_ENCRYPTION_PHASED_VCF } from '../modules/local/compression/compression_encryption_phased_vcf'

workflow ENCRYPTION {
    take:
    chunks

    main:
    if (params.mode == 'phasing_only') {
            COMPRESSION_ENCRYPTION_PHASED_VCF(chunks)
        } else {
            COMPRESSION_ENCRYPTION_VCF(chunks)
        }
}
