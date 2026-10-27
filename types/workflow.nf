record VcfMetafile {
    chr: String
    start: String
    end: String
    phasing_status: String
    vcf: Path
}

record ImputedChunk {
    chr: String
    start: String
    end: String
    dose: Path
    info: Path
    empirical_dose: Path
}

record EagleInput {
    chr: String
    bcf: Path
    bcf_csi: Path
    start: String
    end: String
    phasing_status: String
    vcf: Path
}

record BeagleInput {
    chr: String
    bcf: Path
    start: String
    end: String
    phasing_status: String
    vcf: Path
    map_beagle: Path
}

record MinimacInput {
    chr: String
    start: String
    end: String
    phasing_status: String
    vcf: Path
    m3vcf: Path
}
