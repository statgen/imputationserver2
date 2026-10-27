nextflow.enable.types = true

include { PgsParams } from '../../../types/params'
include { ImputedChunk } from '../../../types/workflow'

process CALCULATE_CHUNKS {
    label 'pgs'
    tag "${chunk.dose}"

    input:
    chunk: ImputedChunk
    scores_txt: Path
    scores_info: Path
    scores_index: Path
    scores: Path
    pgs: PgsParams

    output:
    scores: Path = file("${name}.scores.txt")
    info: Path = file("${name}.scores.info")

    script:
    name = "${chunk.dose.baseName}_${chunk.chr}_${chunk.start}_${chunk.end}"

    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[CALCULATE_CHUNKS] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    """
    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        apply ${chunk.dose} \
        --ref ${scores_txt} \
        --out ${name}.scores.txt \
        --include-scores ${scores} \
        --info ${name}.scores.info \
        --start ${chunk.start} \
        --end ${chunk.end} \
        ${pgs.fix_strand_flips ? "--fix-strand-flips" : ""} \
        --min-r2 ${pgs.min_r2} \
        --no-ansi
    """
}
