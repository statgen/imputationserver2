record PhasingParams {
    engine: String // eagle | beagle
    window: Integer
}

record QcFilterParams {
    overlap: Float?
    minSnps: Integer?
    sampleCallrate: Float?
    mixedGenotypeschrX: Float?
    strandFlips: Integer?
    alleleSwitches: Integer?
}

record PopulationParams {
    id: String
    name: String
    samples: Integer
}

record RefPanelParams {
    id: String
    build: String
    genotypes: String
    sites: String
    refEagle: String?
    mapEagle: String?
    refBeagle: String?
    mapBeagle: String?
    mapMinimac: String?
    qcFilter: QcFilterParams?
    populations: List<PopulationParams>
    range: String?
}

record EncryptionParams {
    enabled: Boolean
    aes: Boolean
    annotate: Boolean
    thread_scale: Integer
}

record ImputationParams {
    enabled: Boolean
    window: Integer
    minimac_min_ratio: Float
    min_r2: Float
    decay: Float
    prob_threshold: Float
    prob_threshold_s1: Float
    diff_threshold: Float
    min_recom: Float
    meta: Boolean
    md5: Boolean
    create_index: Boolean
}

record ServiceParams {
    name: String
    contact: String?
    email: String?
    url: String?
    threads: Integer
}

record UserParams {
    name: String?
    email: String?
}

record AncestryParams {
    enabled: Boolean
    max_pcs: Integer
    k: Integer
    threshold: Float
    dim: Integer
    dim_high: Integer
    batch_size: Integer
    references: String?
}

record PgsParams {
    enabled: Boolean
    category: String
    fix_strand_flips: Boolean
    min_r2: Float
}

record PgsCatalogParams {
    scores: String
    meta: String
}
