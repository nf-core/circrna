include { TARGETSCAN_DATABASE } from '../../modules/local/targetscan/database/main'
include { TARGETSCAN          } from '../../modules/local/targetscan/predict/main'
include { MIRANDA             } from '../../modules/nf-core/miranda/main'
include { MIRNA_TARGETS       } from '../../modules/local/mirna_targets/main'

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

    consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
    MIRNA_TARGETS( consolidate_targets )

    ch_versions = ch_versions.mix(MIRNA_TARGETS.out.versions)

    emit:
    versions = ch_versions
}
