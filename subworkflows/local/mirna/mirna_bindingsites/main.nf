include { BIOAWK as ADD_BACKSPLICE        } from '../../../../modules/nf-core/bioawk'
include { BIOAWK as FORMAT_TARGETSCAN     } from '../../../../modules/nf-core/bioawk'
include { MIRANDA                         } from '../../../../modules/nf-core/miranda'
include { GAWK as UNIFY_MIRANDA           } from '../../../../modules/nf-core/gawk'
include { TARGETSCAN                      } from '../../../../modules/local/targetscan/predict'
include { GAWK as UNIFY_TARGETSCAN        } from '../../../../modules/nf-core/gawk'
include { MIRNA_TARGETS                   } from '../../../../modules/local/mirna_targets'
include { MAJORITY_VOTE                   } from '../../../../modules/local/majority_vote'

workflow MIRNA_BINDINGSITES {
    take:
    transcriptome_fasta
    _circrna_bed12
    mirna_fasta
    tools
    min_tools

    main:
    ch_versions = channel.empty()
    ch_predictions = channel.empty()

    // miRNAs can potentially bind to circRNAs right at the backsplice site
    // In this case, the miRNA binding sequence would partially overlap with start and end of the circRNA
    // To account for this, the first 25bp of the circRNA are added to the end of the circRNA sequence
    ADD_BACKSPLICE( transcriptome_fasta, "fa" )

    ch_transcriptome_batches = ADD_BACKSPLICE.out.output
        .splitFasta(by: 100, file: true)
        .map{ _meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file]}

    //
    // MIRNA PREDICTION TOOLS:
    //
    if (tools.size() == 0) {
        error 'No tools selected for miRNA discovery.'
    }

    if (tools.contains('targetscan')) {
        //
        // TARGETSCAN WORKFLOW:
        //
        // TargetScan expects one line per miRNA: id, seed sequence (7bp after the first base) and a dummy species id
        FORMAT_TARGETSCAN( mirna_fasta, "txt" )
        TARGETSCAN( ch_transcriptome_batches, FORMAT_TARGETSCAN.out.output.collect() )
        UNIFY_TARGETSCAN( TARGETSCAN.out.txt.map { meta, txt -> [meta, txt, 'targetscan.tsv'] }, [], false )

        ch_versions = ch_versions.mix(TARGETSCAN.out.versions)
        ch_predictions = ch_predictions.mix(UNIFY_TARGETSCAN.out.output)
    }

    if (tools.contains('miranda')) {
        //
        // MIRANDA WORKFLOW:
        //
        MIRANDA( ch_transcriptome_batches, mirna_fasta.map{_meta, mature -> mature}.collect() )
        UNIFY_MIRANDA( MIRANDA.out.txt.map { meta, txt -> [meta, txt, 'miranda.tsv'] }, [], false )

        ch_versions = ch_versions.mix(MIRANDA.out.versions)
        ch_predictions = ch_predictions.mix(UNIFY_MIRANDA.out.output)
    }

    //
    // CONSOLIDATE PREDICTIONS WORKFLOW:
    //
    // TODO: This is an artifact and should be removed if we have a replacement

    if (tools.contains('targetscan') && tools.contains('miranda')) {
        // consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt).join(circrna_bed12)
        consolidate_targets = TARGETSCAN.out.txt.join(MIRANDA.out.txt)

        MIRNA_TARGETS( consolidate_targets )

        ch_versions = ch_versions.mix(MIRNA_TARGETS.out.versions)
    }

    //
    // MAJORITY VOTING:
    //
    MAJORITY_VOTE( ch_predictions.map{ _meta, file -> file }.collect().map{ files -> [[id: "mirna"], files] }, min_tools )
    ch_versions = ch_versions.mix(MAJORITY_VOTE.out.versions)

    emit:
    binding_sites = MAJORITY_VOTE.out.targets

    versions = ch_versions
}
