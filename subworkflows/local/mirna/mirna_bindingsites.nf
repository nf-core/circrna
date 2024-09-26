include { BIOAWK as ADD_BACKSPLICE        } from '../../../modules/nf-core/bioawk'
include { CAT_CAT as COMBINE_BINDINGSITES } from '../../../modules/nf-core/cat/cat'
include { GAWK as UNIFY_MIRANDA           } from '../../../modules/nf-core/gawk'
include { GAWK as UNIFY_TARGETSCAN        } from '../../../modules/nf-core/gawk'
include { GAWK as UNIFY_TARPMIR           } from '../../../modules/nf-core/gawk'
include { MAJORITY_VOTE                   } from '../../../modules/local/majority_vote'
include { MIRANDA                         } from '../../../modules/nf-core/miranda'
include { MIRNA_TARGETS                   } from '../../../modules/local/mirna_targets'
include { PITA                            } from '../../../modules/local/pita'
include { TARGETSCAN                      } from '../../../modules/local/targetscan/predict'
include { TARPMIR                         } from '../../../modules/local/tarpmir'

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
    ADD_BACKSPLICE( transcriptome_fasta )
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    ch_transcriptome_batches = ADD_BACKSPLICE.out.output
        .splitFasta(by: 100, file: true)
        .map{ meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file]}

    //
    // MIRNA PREDICTION TOOLS:
    //
    tools_selected = params.mirna_tools.split(',').collect{it.trim().toLowerCase()}

    if (tools_selected.size() == 0) {
        error 'No tools selected for miRNA discovery.'
    }

    if (tools_selected.contains('targetscan')) {
        //
        // TARGETSCAN WORKFLOW:
        //
        TARGETSCAN( ch_transcriptome_batches, formatMiRNAForTargetScan( mirna_fasta ).collect() )
        UNIFY_TARGETSCAN( TARGETSCAN.out.txt, [] )

        ch_versions = ch_versions.mix(TARGETSCAN.out.versions)
        ch_versions = ch_versions.mix(UNIFY_TARGETSCAN.out.versions)
        ch_predictions = ch_predictions.mix(UNIFY_TARGETSCAN.out.output)
    }

    if (tools_selected.contains('miranda')) {
        //
        // MIRANDA WORKFLOW:
        //
        MIRANDA( ch_transcriptome_batches, mirna_fasta.map{meta, mature -> mature}.collect() )
        UNIFY_MIRANDA( MIRANDA.out.txt, [] )

        ch_versions = ch_versions.mix(MIRANDA.out.versions)
        ch_versions = ch_versions.mix(UNIFY_MIRANDA.out.versions)
        ch_predictions = ch_predictions.mix(UNIFY_MIRANDA.out.output)
    }

    if (tools_selected.contains('tarpmir')) {
        //
        // TARPMIR WORKFLOW:
        //

        TARPMIR( ch_transcriptome_batches, mirna_fasta.collect() )
        UNIFY_TARPMIR( TARPMIR.out.bindings, [] )

        ch_versions = ch_versions.mix(TARPMIR.out.versions)
        ch_versions = ch_versions.mix(UNIFY_TARPMIR.out.versions)
        ch_predictions = ch_predictions.mix(UNIFY_TARPMIR.out.output)
    }

    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //
    // TODO: This is an artifact and should be removed if we have a replacement

    // consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
    consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt)

    MIRNA_TARGETS( consolidate_targets )

    ch_versions = ch_versions.mix(MIRNA_TARGETS.out.versions)

    //
    // MAJORITY VOTING:
    //
    MAJORITY_VOTE( ch_predictions.map{meta, file -> file}.collect().map{[[id: "mirna"], it]} )
    ch_versions = ch_versions.mix(MAJORITY_VOTE.out.versions)

    emit:
    targets       = MAJORITY_VOTE.out.targets
    binding_sites = MAJORITY_VOTE.out.binding_sites

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
        .splitFasta(record: [id: true, seqString: true])
        .map { record ->
            return "${record.id}\t${record.seqString[1..7]}\t0000\n"
        }
        .collectFile(name: 'mature.txt')

    ch_targetscan_meta_formatted = ch_targetscan_meta_formatted.map { [[id: "mature_targetscan"], it] }
    return ch_targetscan_meta_formatted
}
