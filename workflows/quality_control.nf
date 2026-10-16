nextflow.enable.types = true

include { RefPanelParams ; PhasingParams } from '../types/params'

include { QUALITY_CONTROL_VCF } from '../modules/local/quality_control/quality_control_vcf'
include { QUALITY_CONTROL_REPORT } from '../modules/local/quality_control/quality_control_report'

workflow QUALITY_CONTROL {
    take:
    validated_files: Bag<Path>
    validation_report: Path
    site_files: Bag<Path>
    refpanel: RefPanelParams
    phasing: PhasingParams
    build: String
    population: String
    chunksize: Integer

    main:
    def chain_file: Path? = null
    if (build != refpanel.build) {
        chain_file = file(
            "${projectDir}/files/chains/${build}To${refpanel.build}.over.chain.gz",
            checkIfExists: true,
        )
    }

    qc = QUALITY_CONTROL_VCF(
        validated_files,
        site_files,
        chain_file,
        refpanel,
        phasing,
        population,
        build,
        chunksize,
    )

    chunks_vcf_index = qc.chunks_vcf
        .flatMap { paths -> paths } // Convert from Value<Set<Path>> to Channel<Path>
        .map { path -> record(key: path.baseName, vcf: path) }

    chunks_csv_index = qc.chunks_csv
        .flatMap { paths -> paths } // Convert from Value<Set<Path>> to Channel<Path>
        .flatMap { path -> path.splitCsv(header: false, sep: '\t') as List<List<String>> }
        .map { row -> record(
            key: file(row[4]).baseName,
            chr: row[0],
            start: row[1],
            end: row[2],
            phasing_status: row[3],
            // snps: row[5],
            // in_reference: row[6],
        ) }

    metafiles_ch = chunks_csv_index
        .join(chunks_vcf_index, by: 'key') // JOIN in the SQL sense.
        .map { row -> tuple(row.chr, row.start, row.end, row.phasing_status, row.vcf) }

    QUALITY_CONTROL_REPORT(
        qc.maf_file,
        validation_report,
        qc.qc_report,
        file("${baseDir}/files/qc-report.Rmd", checkIfExists: true),
        population,
    )

    emit:
    qc_metafiles = metafiles_ch
}
