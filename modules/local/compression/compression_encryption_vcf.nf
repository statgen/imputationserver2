process COMPRESSION_ENCRYPTION_VCF {
    label 'postprocessing'
    publishDir params.output, mode: 'copy'
    tag "Merge Chromosome ${chr}"

    input:
    tuple val(chr), val(start), val(end), path(imputed_vcf_data), path(imputed_info), path(imputed_meta_vcf_data)

    output:
    path "*.zip", emit: encrypted_file, optional: true
    path "*.md5", emit: md5_file, optional: true
    path "chr${chr}*", emit: raw_files, optional: true

    script:
    def output_format = params.imputation.output_format ?: 'vcf.gz'
    if (!['sav', 'vcf.gz'].contains(output_format)) {
            throw new IllegalArgumentException("Invalid params.imputation.output_format: ${output_format}. Use 'sav' or 'vcf.gz'.")
        }

    def imputed_joined = processFileList(imputed_vcf_data)
    def meta_joined = processFileList(imputed_meta_vcf_data)
    def info_joined = processFileList(imputed_info)

    def prefix = "chr${chr}"
    def imputed_name = "${prefix}.dose.${output_format}"
    def meta_name = "${prefix}.empiricalDose.${output_format}"
    def zip_name = "chr_${chr}.zip"
    def info_name = "${prefix}.info.gz"

    def aes = params.encryption.aes ? "-mem=AES256" : ""
    def panel_version = params.refpanel.id
    def scale   = params.encryption.thread_scale ?: 1
    def threads = Math.max(1, Math.floor(task.cpus * scale) as int)
    def first_file = getFirstFile(imputed_vcf_data)

    """
    if [[ "${output_format}" == "sav" ]]
    then
        # SAV format info file
        sav concat  ${info_joined} -o ${info_name}

        if [[ "${params.encryption.annotate}" = true ]]
        then
            sav head -O vcf "${first_file}" | \\
            awk -v pipeline="${workflow.manifest.version}" \\
                -v phasing="${params.phasing.engine}" \\
                -v panel="${panel_version}" '
            {
                print
            }
            NR==1 {
                print "##mis_pipeline=" pipeline
                print "##mis_phasing=" phasing
                print "##mis_panel=" panel
                }
                '  > final_header.vcf

             #bcftools concat --threads ${task.cpus} -n ${imputed_joined} -o intermediate_${imputed_name} -Oz
             sav concat  ${imputed_joined} -o intermediate_${imputed_name}
             sav rehead final_header.vcf intermediate_${imputed_name} -o ${imputed_name}
             rm -f intermediate_${imputed_name}
        else
            sav concat ${imputed_joined} -o ${imputed_name}
        fi

        if [[ "${params.imputation.meta}" = true ]]
        then
            sav concat ${meta_joined} -o ${meta_name}
        fi

        if [[ "${params.imputation.create_index}" = true ]]
        then
            sav index ${imputed_name}
        fi

    else
        # concat info files
        bcftools concat --threads ${threads} -n ${info_joined} -o ${info_name} -Oz

        # annotate files
        if [[ "${params.encryption.annotate}" = true ]]
        then
            # Extracting the original header then paste the new one
            bcftools view -h "${first_file}" > header_from_chunk.txt

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
    fi

    # zip files
    if [[ "${params.encryption.enabled}" = true ]]
    then
        7z a -tzip ${aes} -mmt${threads} -p"${params.encryption_password}" ${zip_name} ${prefix}*
        if [[ "${output_format}" == "sav" ]]
        then
            rm -f *.sav *.s1r *.s1r.gz *.s1r.gz.tbi *info.gz *info.gz.tbi
        else
            rm -f *vcf.gz* *info.gz *info.gz.tbi
        fi
    fi

    # create md5 of zip file
    if [[ "${params.encryption.enabled}" = true && "${params.imputation.md5}" = true ]]
    then
        md5sum ${zip_name} > ${zip_name}.md5
    fi

    # create md5 of imputed file
    if [[ "${params.encryption.enabled}" = false && "${params.imputation.md5}" = true ]]
    then
        md5sum ${imputed_name} > ${imputed_name}.md5
    fi

    """
}

def compareFilenames(a, b) {
    a = a.toString().replaceAll('PAR1', '1').replaceAll('nonPAR', '2').replaceAll('PAR2', '3')
    b = b.toString().replaceAll('PAR1', '1').replaceAll('nonPAR', '2').replaceAll('PAR2', '3')
    return a <=> b
}

def processFileList(fileList) {
    return fileList.sort { a, b -> compareFilenames(a, b) }.join(" ")
}

def getFirstFile(fileList) {
    def sorted = fileList.sort { a, b -> compareFilenames(a, b) }
    return sorted[0]
}

