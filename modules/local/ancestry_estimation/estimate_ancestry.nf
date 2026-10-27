nextflow.enable.types = true

include { AncestryParams } from '../../../types/params'

process ESTIMATE_ANCESTRY {
    label 'ancestry'
    publishDir params.output, mode: 'copy'

    input:
    study_pcs: Bag<Path>
    reference_pc_coord: Path
    reference_samples: Path
    ancestry: AncestryParams

    output:
    populations: Path = file("estimated-population.txt")

    script:
    def study_pcs_joined = study_pcs.join(' ')

    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[ESTIMATE_ANCESTRY] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    """
    # merge csv files
    csvtk concat ${study_pcs_joined} > study.ProPC.coord

    java -Xmx${avail_mem}M -jar /opt/imputationserver-utils/imputationserver-utils.jar \
        estimate-ancestry                    \
        --samples ${reference_samples}       \
        --reference-pc ${reference_pc_coord} \
        --study-pc study.ProPC.coord         \
        --max-pcs ${ancestry.max_pcs}        \
        --k ${ancestry.k}                    \
        --threshold ${ancestry.threshold}    \
        --output estimated-population.txt
    """
}
