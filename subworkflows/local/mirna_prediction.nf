include { TARGETSCAN_DATABASE } from '../../modules/local/targetscan/database/main'
include { TARGETSCAN          } from '../../modules/local/targetscan/predict/main'
include { MIRANDA             } from '../../modules/nf-core/miranda/main'

workflow MIRNA_PREDICTION{

    take:
    circrna_fasta
    circrna_bed12
    mature

    main:
    ch_versions = Channel.empty()

    //
    // TARGETSCAN WORKFLOW:
    //

    TARGETSCAN_DATABASE( mature )
    TARGETSCAN( circrna_fasta, TARGETSCAN_DATABASE.out.mature_txt )

    ch_versions = ch_versions.mix(TARGETSCAN.out.versions)

    //
    // MIRANDA WORKFLOW:
    //

    MIRANDA( circrna_fasta, mature )

    ch_versions = ch_versions.mix(MIRANDA.out.versions)

    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //

    test = TARGETSCAN.out.map{ meta, it -> sample_id = it.simpleName; return [ meta, "${sample_id}.targetscan.txt" ]}.view()

    emit:
    foo = "foo"
}
