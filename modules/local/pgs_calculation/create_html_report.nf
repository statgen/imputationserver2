nextflow.enable.types = true

process CREATE_HTML_REPORT {
    label 'pgs'
    publishDir params.output, mode: 'copy'

    input:
    merged_score: Path
    merged_info: Path
    scores_meta: Path
    estimated_ancestry: Path?

    output:
    html_report: Path = file("scores.html")
    zip_report: Path = file("scores.report.zip")
    coverage_report: Path = file("scores.coverage")

    script:
    def samples = (estimated_ancestry != null) ? "--samples ${estimated_ancestry}" : ""

    def avail_mem: Integer = 1024
    if (!task.memory) {
        log.info('[CREATE_HTML_REPORT] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.toMega() * 0.8).intValue()
    }

    """
    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        report \
        --info ${merged_info} \
        --meta ${scores_meta} \
        ${samples} \
        --out scores.html

    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        report \
        --data ${merged_score} \
        --info ${merged_info} \
        --meta ${scores_meta} \
        --template multi-page \
        ${samples} \
        --out scores.details.html

    7z a scores.report.zip scores.details.html scores
    rm -rf scores.details.html scores

    java -Xmx${avail_mem}M -jar /opt/pgs-calc/pgs-calc.jar \
        report \
        --data ${merged_score} \
        --info ${merged_info} \
        --meta ${scores_meta} \
        --template txt \
        --out scores.coverage
    """
}
