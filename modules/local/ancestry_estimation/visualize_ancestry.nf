nextflow.enable.types = true

process VISUALIZE_ANCESTRY {
    label 'ancestry'
    publishDir params.output, mode: 'copy'

    input:
    ancestry_estimation_report: Path
    estimated_ancestry: Path
    reference_pc_coord: Path
    reference_samples: Path

    output:
    visualization: Set<Path> = files("*.html")

    script:
    """
    # create PCA plot
    Rscript -e "require( 'rmarkdown' ); render('${ancestry_estimation_report}',
        params = list(
            populations = '${estimated_ancestry}',
            reference_pc_coord = '${reference_pc_coord}',
            reference_samples = '${reference_samples}'
        ),
        intermediates_dir='\$PWD',
        knit_root_dir='\$PWD',
        output_file='\$PWD/estimated-populations.html'
    )"
    """
}
