nextflow.enable.types =  true

include { AncestryParams } from '../types/params'

include { PREPARE_TRACE } from '../modules/local/ancestry_estimation/prepare_trace'
include { EXECUTE_TRACE } from '../modules/local/ancestry_estimation/execute_trace'
include { ESTIMATE_ANCESTRY } from '../modules/local/ancestry_estimation/estimate_ancestry'
include { VISUALIZE_ANCESTRY } from '../modules/local/ancestry_estimation/visualize_ancestry'

workflow ANCESTRY_ESTIMATION {
    take:
    vcf_files: Channel<Path> // TODO(Marc): This was renamed from 'files'. Ensure the `vcf_` prefix is accurate.
    ancestry: AncestryParams
    build: String

    main:
    references = files(ancestry.references, glob: true)

    ref_site    = references.find { ref -> ref.getExtension() == 'site'    }
    ref_range   = references.find { ref -> ref.getExtension() == 'range'   }
    ref_geno    = references.find { ref -> ref.getExtension() == 'geno'    }
    ref_coord   = references.find { ref -> ref.getExtension() == 'coord'   }
    ref_samples = references.find { ref -> ref.getExtension() == 'samples' }

    trace = PREPARE_TRACE(vcf_files.collect(), ref_site, ancestry.batch_size, build)

    // Cast Value to Channel (so EXECUTE_TRACE distributes accross paths).
    batch_files = trace.batches.flatMap { paths -> paths }

    pcs = EXECUTE_TRACE(
        batch_files,
        trace.vcf,
        ref_site,
        ref_range,
        ref_geno,
        ref_coord,
        ref_samples,
    )

    populations = ESTIMATE_ANCESTRY(
        pcs,
        ref_coord,
        ref_samples,
        ancestry
    )

    VISUALIZE_ANCESTRY(
        file("${baseDir}/files/ancestry-estimation.Rmd"),
        populations,
        ref_coord,
        ref_samples,
    )

    emit:
    estimated_ancestry = populations
}
