nextflow.enable.types = true

include { RefPanelParams ; ImputationParams ; EncryptionParams } from '../../../types/params'
include { ImputedChunk } from '../../../types/workflow'

process COMPRESSION_ENCRYPTION_VCF {
    label 'postprocessing'
    publishDir params.output, mode: 'copy'
    tag "Merge Chromosome ${chr}"

    input:
    tuple(chr: String, chunks: Bag<ImputedChunk>)
    password: String
    phasing_engine: String
    refpanel: RefPanelParams
    imputation: ImputationParams
    encryption: EncryptionParams

    output:
    encrypted_file: Path? = file("*.zip", optional: true)
    md5_file: Path? = file("*.md5", optional: true)
    raw_files: Set<Path> = files("chr${chr}*", optional: true)

    script:
    def sorted_chunks: List<ImputedChunk> = chunks.toSorted { chunk ->
        chunk.dose.name
            .replaceAll('PAR1', '1')
            .replaceAll('nonPAR', '2')
            .replaceAll('PAR2', '3')
    }

    def imputed_joined = sorted_chunks.collect { chunk -> chunk.dose }.join(' ')
    def info_joined = sorted_chunks.collect { chunk -> chunk.info }.join(' ')
    def meta_joined = sorted_chunks.collect { chunk -> chunk.empirical_dose }.join(' ')

    def prefix = "chr${chr}"
    def imputed_name = "${prefix}.dose.vcf.gz"
    def meta_name = "${prefix}.empiricalDose.vcf.gz"
    def zip_name = "chr_${chr}.zip"
    def info_name = "${prefix}.info.gz"
    def aes = encryption.aes ? '-mem=AES256' : ''
    def panel_version = refpanel.id
    def scale = encryption.thread_scale ?: 1
    def threads = Math.max(1, task.cpus * scale)
    def first_vcf = sorted_chunks[0].dose

    """
    # concat info files
    bcftools concat --threads ${threads} -n ${info_joined} -o ${info_name} -Oz

    # annotate files
    if [[ "${encryption.annotate}" = true ]]
    then
        # Extracting the original header then paste the new one
        bcftools view -h "${first_vcf}" > header_from_chunk.txt

        # build final header: keep ## lines, add custom lines, then #CHROM last
        awk -v pipeline="${workflow.manifest.version}" \
                -v phasing="${params.phasing.engine}" \
                -v panel="${panel_version}" '
            /^#CHROM/ {
                print "##mis_pipeline=" pipeline
                print "##mis_phasing=" phasing
                print "##mis_panel=" panel
            }
            { print }
            ' header_from_chunk.txt > final_header.txt

        bcftools concat --threads ${threads} -n ${imputed_joined} -Oz | bcftools reheader -h final_header.txt -o ${imputed_name}
    else
        bcftools concat --threads ${threads} -n ${imputed_joined} -o intermediate_${imputed_name} -Oz
        mv intermediate_${imputed_name} ${imputed_name}
    fi

    # write meta files
    if [[ "${params.imputation.meta}" = true ]]
    then
        bcftools concat --threads ${threads} -n ${meta_joined} -o ${meta_name} -Oz
        tabix ${meta_name}
    fi

    # create tabix files
    if [[ "${params.imputation.create_index}" = true ]]
    then
        tabix ${imputed_name}
    fi

    # zip files
    if [[ "${encryption.enabled}" = true ]]
    then
        7z a -tzip ${aes} -mmt${threads} -p"${password}" ${zip_name} ${prefix}*
        rm *vcf.gz* *info.gz
    fi

    # create md5 of zip file
    if [[ "${encryption.enabled}" = true && "${params.imputation.md5}" = true ]]
    then
        md5sum ${zip_name} > ${zip_name}.md5
    fi

    # create md5 of imputed file
    if [[ "${encryption.enabled}" = false && "${params.imputation.md5}" = true ]]
    then
        md5sum ${imputed_name} > ${imputed_name}.md5
    fi

    """
}
