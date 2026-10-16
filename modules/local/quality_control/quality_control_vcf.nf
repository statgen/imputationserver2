nextflow.enable.types = true

include { RefPanelParams ; PhasingParams } from '../../../types/params'

process QUALITY_CONTROL_VCF {
    label 'preprocessing'
    publishDir params.output, mode: 'copy', pattern: "qc_report.txt"
    publishDir params.output, mode: 'copy', pattern: "statistics/*.txt"

    input:
    vcf_files: Bag<Path>
    site_files: Bag<Path> // Not used explicitly, but adding as input causes files to be staged.
    chain_file: Path?
    refpanel: RefPanelParams
    phasing: PhasingParams
    population: String
    build: String
    chunksize: Integer

    output:
    chunks_csv: Set<Path> = files("metafiles/*", optional: true)
    chunks_vcf: Set<Path> = files("chunks/*")
    statistics: Set<Path> = files("statistics/*", optional: true)
    maf_file: Path? = file("maf.txt", optional: true)
    qc_report: Path = file("qc_report.txt")

    script:
    def avail_mem = 1024
    if (!task.memory) {
        log.info('[QUALITY_CONTROL_VCF] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    def chain = (chain_file != null) ? "--chain ${chain_file}" : ''

    def refpanel_json = groovy.json.JsonOutput.toJson(refpanel)

    // Typed processes interpolate collections as `[a, b, c]` with brackets (new behavior).
    // Hence, we need to explicitly turn vcf_files into a space separated list.
    def vcf_args = vcf_files.join(' ')

    """
    set +e
    echo '${refpanel_json}' > reference-panel.json

    # Verify if VCF files are valid
    for vcf in ${vcf_args}; do
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

    java -Xmx${avail_mem}M -jar /opt/imputationserver-utils/imputationserver-utils.jar \
        run-qc \
        --reference reference-panel.json \
        --population ${population} \
        --build ${build} \
        --chunksize ${chunksize} \
        --phasing-window ${phasing.window} \
        --maf-output maf.txt \
        --chunks-out chunks \
        --statistics-out statistics \
        --metafiles-out metafiles \
        --report qc_report.txt \
        --no-index \
        ${chain} \
        ${vcf_args}

    exit_code_a=\$?

    # Check if QC step failed
    if [[ \$exit_code_a -ne 0 ]]; then
        rm -rf metafiles
    fi

    cat qc_report.txt

    # Always succeed, so QC files get published.
    exit 0
    """
}
