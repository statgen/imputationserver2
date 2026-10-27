nextflow.enable.types = true

include { RefPanelParams ; PhasingParams } from '../types/params'
include { VcfMetafile } from '../types/workflow'

include { EAGLE } from '../modules/local/phasing/eagle'
include { BEAGLE } from '../modules/local/phasing/beagle'

workflow PHASING {
    take:
    metafiles: Channel<VcfMetafile>
    refpanel: RefPanelParams
    phasing: PhasingParams
    threads: Integer

    main:
    def phased_ch: Channel<VcfMetafile>

    def chromosomes = channel.fromList(1..22)
        .map { chr -> "${chr}" }
        .mix(channel.of('X.nonPAR', 'X.PAR1', 'X.PAR2', 'MT'))

    if (phasing.engine == 'eagle') {
        if (refpanel.refEagle == null || refpanel.mapEagle == null) {
            error('Eagle requires refpanel.refEagle and refpanel.mapEagle')
        }

        phasing_map = file(refpanel.mapEagle, checkIfExists: true)

        def eagle_references = chromosomes
            .map { chr ->
                record(
                    chr: chr,
                    bcf: file(PatternUtil.parse(refpanel.refEagle, [chr: chr])),
                    bcf_csi: file(PatternUtil.parse(refpanel.refEagle + ".csi", [chr: chr])),
                )
            }
            .filter { ref -> ref.bcf.exists() && ref.bcf_csi.exists() }

        def eagle_bcf_metafiles_ch = eagle_references.join(metafiles, by: 'chr')

        phased_ch = EAGLE(
            eagle_bcf_metafiles_ch,
            phasing_map,
            refpanel.build,
            phasing.window,
            threads,
        )
    }
    else if (phasing.engine == 'beagle') {
        if (refpanel.refBeagle == null || refpanel.mapBeagle == null) {
            error('Beagle requires refpanel.refBeagle and refpanel.mapBeagle')
        }

        def beagle_references = chromosomes
            .map { chr ->
                record(
                    chr: chr,
                    bcf: file(PatternUtil.parse(refpanel.refBeagle, [chr: chr])),
                    map_beagle: file(PatternUtil.parse(refpanel.mapBeagle, [chr: chr])),
                )
            }
            .filter { ref -> ref.bcf.exists() && ref.map_beagle.exists() }

        def beagle_bcf_metafiles_ch = beagle_references.join(metafiles, by: 'chr')

        phased_ch = BEAGLE(
            beagle_bcf_metafiles_ch,
            refpanel.build,
        )
    }
    else {
        error("Unrecognized phasing engine: ${phasing.engine}")
    }

    emit:
    phased: Channel<VcfMetafile> = phased_ch
}
