process COMPRESSION_ENCRYPTION_PHASED_VCF {
    label 'postprocessing'
    publishDir params.output, mode: 'copy'
    tag "Merge Phased Chromosome ${chr}"

    input:
    tuple val(chr), val(start), val(end), val(phasing_status), path(phased_vcf_data)

    output:
    path "*.zip", emit: encrypted_file, optional: true
    path "*.md5", emit: md5_file, optional: true
    path "chr${chr}*", emit: raw_files, optional: true

    script:
    def prefix = "chr${chr}"
    def phased_name = "${prefix}.phased.vcf.gz"
    def zip_name = "chr_${chr}.zip"
    def aes = params.encryption.aes ? "-mem=AES256" : ""
    def panel_version = params.refpanel.id
    def scale   = params.encryption.thread_scale ?: 1
    def threads = Math.max(1, Math.floor(task.cpus * scale) as int)
    def entries = [start, end, phased_vcf_data].transpose().sort { a, b ->
        def startCompare = a[0].toString().toLong() <=> b[0].toString().toLong()
        startCompare != 0 ? startCompare : compareFilenames(a[2], b[2])
    }
    def index_commands = entries.collect { entry ->
        def chunk_file = entry[2]
        """
        tabix -f ${chunk_file}
        """
    }.join("\n")
    def phased_joined = entries.collect { it[2] }.join(" ")
    def first_vcf = entries[0][2]

    """
    ${index_commands}

    # annotate files
    if [[ "${params.encryption.annotate}" = true ]]
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

        bcftools concat --threads ${threads} --ligate -Oz ${phased_joined} | bcftools reheader -h final_header.txt -o ${phased_name}
    else
        bcftools concat --threads ${threads} --ligate -Oz -o intermediate_${phased_name} ${phased_joined}
        mv intermediate_${phased_name} ${phased_name}
    fi

    # create tabix files
    if [[ "${params.imputation.create_index}" = true ]]
    then
        tabix ${phased_name}
    fi

    # zip files
    if [[ "${params.encryption.enabled}" = true ]]
    then
        7z a -tzip ${aes} -mmt${threads} -p"${params.encryption_password}" ${zip_name} ${prefix}*
        rm *vcf.gz*
    fi

    # create md5 of zip file
    if [[ "${params.encryption.enabled}" = true && "${params.imputation.md5}" = true ]]
    then
        md5sum ${zip_name} > ${zip_name}.md5
    fi

    # create md5 of phased file
    if [[ "${params.encryption.enabled}" = false && "${params.imputation.md5}" = true ]]
    then
        md5sum ${phased_name} > ${phased_name}.md5
    fi

    """
}

def compareFilenames(a, b) {
    a = a.toString().replaceAll('PAR1', '1').replaceAll('nonPAR', '2').replaceAll('PAR2', '3')
    b = b.toString().replaceAll('PAR1', '1').replaceAll('nonPAR', '2').replaceAll('PAR2', '3')
    return a <=> b
}
