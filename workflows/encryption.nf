nextflow.enable.types = true

include { RefPanelParams ; ImputationParams ; EncryptionParams } from '../types/params'
include { ImputedChunk } from '../types/workflow'

include { COMPRESSION_ENCRYPTION_VCF } from '../modules/local/compression/compression_encryption_vcf'

workflow ENCRYPTION {
    take:
    imputed_chunks: Channel<ImputedChunk>
    password: String
    phasing_engine: String
    refpanel: RefPanelParams
    imputation: ImputationParams
    encryption: EncryptionParams

    main:
    def grouped_chunks = imputed_chunks
        .map { chunk -> tuple(chunk.chr, chunk) } // groupBy takes pairs: (key, value)
        .groupBy()

    COMPRESSION_ENCRYPTION_VCF(grouped_chunks, password, phasing_engine, refpanel, imputation, encryption)
}
