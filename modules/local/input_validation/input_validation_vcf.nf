nextflow.enable.types = true

include { RefPanelParams ; ImputationParams ; PhasingParams } from '../../../types/params'

process INPUT_VALIDATION_VCF {
    label 'preprocessing'
    publishDir params.output, mode: 'copy', pattern: '*.{html,log}'

    input:
    population: String
    vcf_files: Bag<Path>
    refpanel: RefPanelParams
    imputation: ImputationParams
    phasing: PhasingParams
    build: String
    mode: String
    chunksize: Integer
    min_samples: Integer
    max_samples: Integer

    output:
    validated_files: Set<Path> = files("*.vcf.gz", includeInputs: true)
    validation_report: Path = file("validation_report.txt")

    script:
    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[INPUT_VALIDATION_VCF] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

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

    java -Xmx${avail_mem}M -jar /opt/imputationserver-utils/imputationserver-utils.jar \
        validate \
        --reference reference-panel.json \
        --population ${population} \
        --build ${build} \
        --r2Filter ${imputation.min_r2} \
        --phasing ${phasing.engine} \
        --mode ${mode} \
        --chunksize ${chunksize} \
        --minSamples ${min_samples} \
        --maxSamples ${max_samples} \
        --report validation_report.txt \
        --no-index \
        ${vcf_args}

    exit_code_a=\$?

    cat validation_report.txt

    exit \$exit_code_a
    """
}
