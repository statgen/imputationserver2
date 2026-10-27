nextflow.enable.types = true

include { RefPanelParams ; PhasingParams } from '../../../types/params'
include { VcfMetafile ; EagleInput } from '../../../types/workflow'

process EAGLE {
    label 'phasing'
    tag "${in.vcf}"

    input:
    in: EagleInput
    map_eagle: Path
    build: String
    phasing_window: Integer
    threads: Integer

    output:
    phased: VcfMetafile = record(
        chr: in.chr,
        start: in.start,
        end: in.end,
        phasing_status: in.phasing_status,
        vcf: file("*.phased.vcf.gz"),
    )

    script:
    // Essentially baseName (but we remove two extensions and we accept names with dots: X.*)
    def chunkfile_name = in.vcf.name.replaceAll('.vcf.gz', '')
    // replace X.nonPAR etc with X for phasing
    def chr_clean = in.chr.startsWith('X.') ? 'X' : in.chr
    def chr_mapped = build == 'hg38' ? 'chr' + chr_clean : chr_clean
    def phasing_start = in.start.toLong() - phasing_window
    phasing_start = phasing_start < 0 ? 1 : phasing_start
    def phasing_end = in.end.toLong() + phasing_window
    def used_threads = threads != -1 ? threads : task.cpus

    """
    tabix ${in.vcf}
    eagle \
        --vcfRef ${in.bcf}  \
        --vcfTarget ${in.vcf} \
        --geneticMapFile ${map_eagle} \
        --outPrefix ${chunkfile_name}.phased \
        --chrom ${chr_mapped} \
        --bpStart ${phasing_start} \
        --bpEnd ${phasing_end} \
        --allowRefAltSwap \
        --vcfOutFormat z \
        --keepMissingPloidyX \
        --numThreads ${used_threads}
    """
}
