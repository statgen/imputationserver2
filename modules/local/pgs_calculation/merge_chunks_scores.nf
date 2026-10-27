nextflow.enable.types = true

process MERGE_CHUNKS_SCORES {
    label 'pgs'
    publishDir params.output, mode: 'copy'

    input:
    score_chunks: Bag<Path>

    output:
    merged_score_files: Path = file("*.txt")

    script:
    def score_chunks_joined = score_chunks.toSorted { path -> path.name }.join(' ')

    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[MERGE_CHUNKS_SCORES] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    """
    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        merge-score ${score_chunks_joined} \
        --out scores.txt
    """
}
