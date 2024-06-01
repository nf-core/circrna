include { TARGETSCAN_DATABASE             } from '../../modules/local/targetscan/database'
include { TARGETSCAN                      } from '../../modules/local/targetscan/predict'
include { MIRANDA                         } from '../../modules/nf-core/miranda'
include { MIRNA_TARGETS                   } from '../../modules/local/mirna_targets'
include { DESEQ2_NORMALIZATION            } from '../../modules/local/deseq2/normalization'
include { MIRNA_FILTERING                 } from '../../modules/local/mirna_filtering'
include { GAWK as MIRANDA_TO_MAJORITY     } from '../../modules/nf-core/gawk'
include { GAWK as TARGETSCAN_TO_MAJORITY  } from '../../modules/nf-core/gawk'
include { CAT_CAT as COMBINE_BINDINGSITES } from '../../modules/nf-core/cat/cat'
include { MAJORITY_VOTE                   } from '../../modules/local/majority_vote' 
include { COMPUTE_CORRELATIONS            } from '../../modules/local/compute_correlations' 

workflow MIRNA_PREDICTION{

    take:
    circrna_fasta
    circrna_bed12
    ch_mature
    ch_mirna
    circrna_counts

    main:
    ch_versions = Channel.empty()
    ch_predictions = Channel.empty()
    ch_votings = Channel.empty()
    ch_targetscan_meta_formatted = Channel.empty()

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

    ch_targetscan_meta_formatted = ch_targetscan_meta_formatted.mix(formatMiRNAForTargetScan( ch_mature )).collect()

    TARGETSCAN( circrna_fasta, ch_targetscan_meta_formatted )
    TARGETSCAN_TO_MAJORITY( TARGETSCAN.out.txt, [] )

    ch_versions = ch_versions.mix(TARGETSCAN.out.versions)
    ch_predictions = ch_predictions.mix(TARGETSCAN_TO_MAJORITY.out.output)

    //
    // MIRANDA WORKFLOW:
    //

    MIRANDA( circrna_fasta, ch_mature.map{meta, mature -> mature} )
    MIRANDA_TO_MAJORITY( MIRANDA.out.txt, [] )
    
    ch_versions = ch_versions.mix(MIRANDA.out.versions)
    ch_predictions = ch_predictions.mix(MIRANDA_TO_MAJORITY.out.output)


    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //

    consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
    MIRNA_TARGETS( consolidate_targets )

    ch_versions = ch_versions.mix(MIRNA_TARGETS.out.versions)


    //
    // UNIFICATION OF PREDICTIONS
    //

    ch_combined = ch_predictions.map{meta, file -> file}.collect().map{[[id: "mirna"], it]}
    
    COMBINE_BINDINGSITES ( ch_combined )
    
    ch_votings = ch_votings.mix(COMBINE_BINDINGSITES.out.file_out)
    ch_versions.mix(COMBINE_BINDINGSITES.out.versions)

    //
    // MAJORITY VOTE: 
    //

    ch_bindingsites = Channel.empty()

    MAJORITY_VOTE( ch_votings )

    ch_bindingsites = ch_bindingsites.mix(MAJORITY_VOTE.output.tsv)
    ch_versions = ch_versions.mix(MAJORITY_VOTE.output.versions)

    //
    // COMPUTE CORREALTION:
    //

    COMPUTE_CORRELATIONS( ch_bindingsites, 
                          ch_mirna_filtered.map{meta, file -> file}.collect(),
                          circrna_counts.map{meta, file -> file}.collect()
                        )

    ch_versions = ch_versions.mix(COMPUTE_CORRELATIONS.out.versions) 
    
    emit:
    versions = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
// Fromatting miRNA input for taretscan
// takes mature.fa, iterates over entries (id, seq) and generates a new file
// writing:
// 1. miR ID
// 2. miR (7bp) seed sequence from mature seq
// 3. Species ID (set to 0000, not important for output).
// to new file
def formatMiRNAForTargetScan(ch_mature) {

    def ch_targetscan_meta_formatted = ch_mature
        .map { meta, mature -> mature }
        .splitFasta(by: 1)
        .map { entry ->
            def lines = entry.readLines()
            def id = lines[0]
            def seq = lines[1..-1].join()
            return [id, seq]
        }
        .map { id, seq ->
            def newSeq = seq[1..7]  // Extract sub-sequence (2)
            return [id, newSeq, '0000']  // Add species ID (3)
        }
        .map { id, newSeq, s_id -> "$id\t$newSeq\t$s_id\n" }
        .collectFile(name: 'mature.txt')

    ch_targetscan_meta_formatted = ch_targetscan_meta_formatted.map { [[id: "mature_targetscan"], it] }
    return ch_targetscan_meta_formatted
}
