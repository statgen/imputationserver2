nextflow.enable.types = true

process QUALITY_CONTROL_REPORT {
    label 'preprocessing'
    publishDir params.output, mode: 'copy'

    input:
    maf_file: Path? // TODO: render() might break if population != "off" and maf_file == null.
    validation_report: Path
    qc_report: Path
    qc_report_file: Path
    population: String

    output:
    reports: Set<Path> = files("*.html")

    script:
    """
    Rscript -e "require( 'rmarkdown' ); render('${qc_report_file}',
        params = list(
            maf_file = '${maf_file}',
            validation_report = '${validation_report}',
            qc_report = '${qc_report}',
            name = '${params.project}',
            allele_frequency_population = '${population}',
            version = '${workflow.manifest.version}',
            date = '${params.project_date}',
            service = '${params.service.name}'
        ),
        intermediates_dir='\$PWD',
        knit_root_dir='\$PWD',
        output_file='\$PWD/quality-control.html'
    )"
    """
}
