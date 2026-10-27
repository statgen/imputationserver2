nextflow.enable.types = true

process EXECUTE_TRACE {
    label 'ancestry'

    input:
    samples: Path
    vcf_file: Path
    reference_site: Path
    reference_range: Path
    reference_geno: Path
    reference_pc_coord: Path
    reference_samples: Path
    dim: Integer
    dim_high: Integer

    output:
    pcs: Path = file("${batch_name}.ProPC.coord")

    script:
    batch_name = "batch_${samples.baseName}"

    """
    tabix ${vcf_file}

    # convert to geno
    vcf2geno \
        --inVcf ${vcf_file} \
        --rangeFile ${reference_range} \
        --peopleIncludeFile ${samples} \
        --out ${batch_name}

    # write config file for trace
    echo "GENO_FILE ${reference_geno}" > trace.config
    echo "STUDY_FILE ${batch_name}.geno" >> trace.config
    echo "COORD_FILE ${reference_pc_coord}" >> trace.config
    echo "DIM ${dim}" >> trace.config
    echo "DIM_HIGH ${dim_high}" >> trace.config
    echo "OUT_PREFIX ${batch_name}" >> trace.config
    echo "NUM_THREADS 1" >> trace.config

    # execute trace with config file
    trace -p trace.config
    """
}
