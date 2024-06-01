include { BIOAWK as ADD_BACKSPLICE        } from '../../../modules/nf-core/bioawk'
include { MIRANDA                         } from '../../../modules/nf-core/miranda'
include { GAWK as UNIFY_MIRANDA           } from '../../../modules/nf-core/gawk'
include { TARGETSCAN                      } from '../../../modules/local/targetscan/predict'
include { GAWK as UNIFY_TARGETSCAN        } from '../../../modules/nf-core/gawk'
include { MIRNA_TARGETS                   } from '../../../modules/local/mirna_targets'
include { CAT_CAT as COMBINE_BINDINGSITES } from '../../../modules/nf-core/cat/cat'
include { MAJORITY_VOTE                   } from '../../../modules/local/majority_vote'

workflow MIRNA_BINDINGSITES {
    take:
    transcriptome_fasta
    circrna_bed12
    mirna_fasta

    main:
    ch_versions = Channel.empty()
    ch_predictions = Channel.empty()

    // miRNAs can potentially bind to circRNAs right at the backsplice site
    // In this case, the miRNA binding sequence would partially overlap with start and end of the circRNA
    // To account for this, the first 25bp of the circRNA are added to the end of the circRNA sequence
    ADD_BACKSPLICE( transcriptome_fasta)
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    ch_transcriptome_batches = ADD_BACKSPLICE.out.output
        .splitFasta(by: 100, file: true)
        .map{ meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file]}

    //
    // TARGETSCAN WORKFLOW:
    //
    TARGETSCAN( ch_transcriptome_batches, formatMiRNAForTargetScan( mirna_fasta ).collect() )
    UNIFY_TARGETSCAN( TARGETSCAN.out.txt, [] )

    ch_versions = ch_versions.mix(TARGETSCAN.out.versions)
    ch_versions = ch_versions.mix(UNIFY_TARGETSCAN.out.versions)
    ch_predictions = ch_predictions.mix(UNIFY_TARGETSCAN.out.output)

    //
    // MIRANDA WORKFLOW:
    //

    MIRANDA( ch_transcriptome_batches, mirna_fasta.map{meta, mature -> mature} )
    UNIFY_MIRANDA( MIRANDA.out.txt, [] )

    ch_versions = ch_versions.mix(MIRANDA.out.versions)
    ch_versions = ch_versions.mix(UNIFY_MIRANDA.out.versions)
    ch_predictions = ch_predictions.mix(UNIFY_MIRANDA.out.output)

    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //
    // TODO: This is an artifact and should be removed if we have a replacement

    consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
    MIRNA_TARGETS( consolidate_targets )

    ch_versions = ch_versions.mix(MIRNA_TARGETS.out.versions)

    //
    // MAJORITY VOTING:
    //
    COMBINE_BINDINGSITES ( ch_predictions.map{meta, file -> file}.collect().map{[[id: "mirna"], it]} )
    MAJORITY_VOTE( COMBINE_BINDINGSITES.out.file_out )

    ch_versions = ch_versions.mix(COMBINE_BINDINGSITES.out.versions)
    ch_versions = ch_versions.mix(MAJORITY_VOTE.out.versions)

    emit:
    binding_sites = MAJORITY_VOTE.out.tsv

    versions = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
// Formatting miRNA input for targetscan
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
