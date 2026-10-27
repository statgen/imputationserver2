nextflow.enable.types = true

include { MinimacInput ; ImputedChunk } from '../../../types/workflow'

process MINIMAC4 {
    label 'imputation'
    tag "${in.vcf}"

    input:
    in: MinimacInput
    minimac_map: Path?
    refpanel_build: String
    minimac_window: Integer
    minimac_min_ratio: Float
    min_r2: Float
    decay: Float
    diffThreshold: Float
    probThreshold: Float
    probThresholdS1: Float
    minRecombination: Float
    threads: Integer

    output:
    imputed_chunks: ImputedChunk = record(
        chr: chr_clean,
        start: in.start,
        end: in.end,
        dose: file("*.dose.vcf.gz"),
        info: file("*.info.gz"),
        empirical_dose: file("*.empiricalDose.vcf.gz"),
    )

    script:
    def map = minimac_map ? '--map ' + "${minimac_map}" : ''
    def r2_filter = min_r2 != 0 ? '--min-r2 ' + "${min_r2}" : ''
    def diff_threshold = diffThreshold != -1 ? '--diff-threshold ' + "${diffThreshold}" : ''
    def prob_threshold = probThreshold != -1 ? '--prob-threshold ' + "${probThreshold}" : ''
    def prob_threshold_s1 = probThresholdS1 != -1 ? '--prob-threshold-s1 ' + "${probThresholdS1}" : ''
    def min_recom = minRecombination != -1 ? '--min-recom ' + "${minRecombination}" : ''
    def chunkfile_name = in.vcf.name.replaceAll('.vcf.gz', '')
    chr_clean = in.chr.startsWith('X.') ? 'X' : in.chr
    def chr_mapped = (refpanel_build == 'hg38') ? 'chr' + chr_clean : chr_clean
    def used_threads = (threads != -1) ? threads : task.cpus

    """
    tabix ${in.vcf}

    minimac4 \
        --region ${chr_mapped}:${in.start}-${in.end} \
        --overlap ${minimac_window} \
        --output ${chunkfile_name}.dose.vcf.gz \
        --output-format vcf.gz \
        --format GT,DS,GP,HDS \
        --min-ratio ${minimac_min_ratio} \
        --all-typed-sites \
        --sites ${chunkfile_name}.info.gz \
        --empirical-output ${chunkfile_name}.empiricalDose.vcf.gz \
        --threads ${used_threads} \
        --decay ${decay} \
        --temp-prefix ./ \
        ${diff_threshold} \
        ${prob_threshold} \
        ${prob_threshold_s1} \
        ${min_recom} \
        ${r2_filter} \
        ${map} \
        ${in.m3vcf} \
        ${in.vcf}
    """
}
