include { TARGETSCAN_DATABASE } from '../../modules/local/targetscan/database/main'
include { TARGETSCAN          } from '../../modules/local/targetscan/predict/main'

workflow MIRNA_PREDICTION{

    take:
    circrna_fasta
    mature

    main:

    TARGETSCAN_DATABASE( mature )
    TARGETSCAN( circrna_fasta, TARGETSCAN_DATABASE.out.mature_txt )


    emit:
    foo = "foo"
}
