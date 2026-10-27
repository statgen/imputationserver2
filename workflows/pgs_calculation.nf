nextflow.enable.types = true

include { PgsParams ; PgsCatalogParams } from '../types/params'
include { ImputedChunk } from '../types/workflow'

include { CALCULATE_CHUNKS } from '../modules/local/pgs_calculation/calculate_chunks'
include { MERGE_CHUNKS_INFOS } from '../modules/local/pgs_calculation/merge_chunks_infos'
include { MERGE_CHUNKS_SCORES } from '../modules/local/pgs_calculation/merge_chunks_scores'
include { CREATE_HTML_REPORT } from '../modules/local/pgs_calculation/create_html_report'
include { FILTER_BY_CATEGORY } from '../modules/local/pgs_calculation/filter_by_category'

workflow PGS_CALCULATION {
    take:
    imputed_chunks: Channel<ImputedChunk>
    estimated_ancestry: Value<Path?>
    pgs: PgsParams
    pgsCatalog: PgsCatalogParams

    main:
    scores_txt = file(pgsCatalog.scores, checkIfExists: true)
    scores_info = file("${pgsCatalog.scores}.info", checkIfExists: true)
    scores_index = file("${pgsCatalog.scores}.tbi", checkIfExists: true)
    scores_meta = file(pgsCatalog.meta, checkIfExists: true)

    def filtered_scores = FILTER_BY_CATEGORY(
        scores_meta,
        pgs.category,
    )

    def calculated_chunks = CALCULATE_CHUNKS(
        imputed_chunks,
        scores_txt,
        scores_info,
        scores_index,
        filtered_scores,
        pgs,
    )

    def merged_scores = MERGE_CHUNKS_SCORES(
        calculated_chunks.scores.collect()
    )

    def merged_infos = MERGE_CHUNKS_INFOS(
        calculated_chunks.info.collect()
    )

    CREATE_HTML_REPORT(
        merged_scores,
        merged_infos,
        scores_meta,
        estimated_ancestry,
    )
}
