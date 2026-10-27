nextflow.enable.types = true

include { RefPanelParams } from '../../../types/params'
include { VcfMetafile ; BeagleInput } from '../../../types/workflow'

process BEAGLE {
    label 'phasing'
    tag "${in.vcf}"

    input:
    in: BeagleInput
    build: String

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
    def phasing_start = in.start.toLong() - params.phasing.window
    phasing_start = phasing_start < 0 ? 1 : phasing_start
    def phasing_end = in.end.toLong() + params.phasing.window
    def used_threads = params.service.threads != -1 ? params.service.threads : task.cpus

    """
    java -jar /usr/bin/beagle.jar \
        ref=${in.bcf} \
        gt=${in.vcf} \
        out=${chunkfile_name}.phased \
        nthreads=${used_threads} \
        chrom=${chr_mapped}:${phasing_start}-${phasing_end} \
        map=${in.map_beagle} \
        impute=false
    """
}
