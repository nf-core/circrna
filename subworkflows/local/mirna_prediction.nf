include { TARGETSCAN_DATABASE            } from '../../modules/local/targetscan/database'
include { TARGETSCAN                     } from '../../modules/local/targetscan/predict'
include { MIRANDA                        } from '../../modules/nf-core/miranda'
include { MIRNA_TARGETS                  } from '../../modules/local/mirna_targets'
include { DESEQ2_NORMALIZATION           } from '../../modules/local/deseq2/normalization'
include { MIRNA_FILTERING                } from '../../modules/local/mirna_filtering'
include { GAWK as MIRANDA_TO_MAJORITY    } from '../../modules/nf-core/gawk'
include { GAWK as TARGETSCAN_TO_MAJORITY } from '../../modules/nf-core/gawk'

workflow MIRNA_PREDICTION{

    take:
    circrna_fasta
    circrna_bed12
    ch_mature
    ch_mirna

    main:
    ch_versions = Channel.empty()

    ADD_BACKSPLICE( circrna_fasta, [] )
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    //
    // MIRNA QUANTIFICATION WORKFLOW:
    //

    ch_mirna_normalized = DESEQ2_NORMALIZATION( ch_mirna ).normalized

    ch_versions = ch_versions.mix(DESEQ2_NORMALIZATION.out.versions)

    ch_mirna_filtered = MIRNA_FILTERING( ch_mirna_normalized, 
	                                     params.mirna_min_sample_percentage,  
										 params.mirna_min_reads
										 ).filtered
	ch_versions = ch_versions.mix(MIRNA_FILTERING.out.versions)
    //
    // TARGETSCAN WORKFLOW:
    //

    TARGETSCAN_DATABASE( ch_mature )
    TARGETSCAN( circrna_fasta, TARGETSCAN_DATABASE.out.mature_txt )
	TARGETSCAN_TO_MAJORITY( TARGETSCAN.out.txt, [] )

    ch_versions = ch_versions.mix(TARGETSCAN_DATABASE.out.versions)
    ch_versions = ch_versions.mix(TARGETSCAN.out.versions)

    //
    // MIRANDA WORKFLOW:
    //

    MIRANDA( circrna_fasta, ch_mature.map{meta, mature -> mature} )
	MIRANDA_TO_MAJORITY( MIRANDA.out.txt, [] )
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
