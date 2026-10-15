process QUALITY_CONTROL_VCF {
    label 'preprocessing'
    publishDir params.output, mode: 'copy', pattern: "qc_report.txt"
    publishDir params.output, mode: 'copy', pattern: "statistics/*.txt"

    input:
    path vcf_files
    path site_files
    path chain_file
    val panel_version
    val refpanel
    val population // String

    output:
    path "metafiles/*", emit: chunks_csv, optional: true
    path "chunks/*", emit: chunks_vcf
    path "statistics/*", optional: true
    path "maf.txt", emit: maf_file, optional: true
    path "qc_report.txt", emit: qc_report

    script:
    def chain = (!panel_version.equals(params.build)) ? "--chain ${chain_file}": ''
    def avail_mem = 1024
    if (!task.memory) {
        log.info('[QUALITY_CONTROL_VCF] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    set +e
    echo '${groovy.json.JsonOutput.toJson(refpanel)}' > reference-panel.json

    # Verify if VCF files are valid
    for vcf in ${vcf_files}; do
        # Attempt to create the index using tabix
        if ! output=\$(tabix -p vcf "\$vcf" 2>&1); then
            echo ::group type=error::
            echo "The provided VCF file is malformed."
            echo "Error: \$output"
            echo ::endgroup::
            exit 1
        fi
    done

    # TODO: create directories in java
    mkdir chunks
    mkdir metafiles
    mkdir statistics
    echo ${chain_file}

    java -Xmx${avail_mem}M -jar /opt/imputationserver-utils/imputationserver-utils.jar \
        run-qc                                    \
        --reference reference-panel.json          \
        --population ${population}                \
        --build ${params.build}                   \
        --chunksize ${params.chunksize}           \
        --phasing-window ${params.phasing.window} \
        --maf-output maf.txt                      \
        --chunks-out chunks                       \
        --statistics-out statistics               \
        --metafiles-out metafiles                 \
        --report qc_report.txt                    \
        --no-index                                \
        ${chain}                                  \
        ${vcf_files}

    exit_code_a=\$?

    # Check if QC step failed
    if [[ \$exit_code_a -ne 0 ]]; then
        rm -rf metafiles
    fi

    cat qc_report.txt

    # Always exit 0 that QC files get published
    exit 0
    """
}
