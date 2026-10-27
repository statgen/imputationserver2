nextflow.enable.types = true

process MERGE_CHUNKS_INFOS {
    label 'pgs'
    publishDir params.output, mode: 'copy'

    input:
    report_chunks: Bag<Path>

    output:
    merged_info_files: Path = file("scores.info")

    script:
    def report_chunks_joined = report_chunks.toSorted { path -> path.name }.join(' ')

    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[MERGE_CHUNKS_INFOS] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    """
    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        merge-info ${report_chunks_joined} \
        --out scores.info
    """
}
