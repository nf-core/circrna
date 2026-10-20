// MODULES
include { DESEQ2_NORMALIZATION } from '../../../modules/local/deseq2/normalization'
include { MIRNA_FILTERING      } from '../../../modules/local/mirna_filtering'
include { COMPUTE_CORRELATIONS } from '../../../modules/local/compute_correlations'

// SUBWORKFLOWS
include { MIRNA_BINDINGSITES   } from '../mirna/mirna_bindingsites'

workflow MIRNA_PREDICTION {
    take:
    transcriptome_fasta
    circrna_annotation
    ch_mature
    mirna_expression
    _transcript_counts
    quantification_rds
    mirna_tools
    mirna_min_tools
    mirna_min_reads
    mirna_min_sample_percentage

    main:
    ch_correlations = channel.empty()

    //
    // MIRNA NORMALIZATION WORKFLOW:
    //

    if (mirna_expression) {

        ch_mirna_normalized = DESEQ2_NORMALIZATION(channel.value(mirna_expression)).normalized


        ch_mirna_filtered = MIRNA_FILTERING(
            ch_mirna_normalized,
            mirna_min_sample_percentage,
            mirna_min_reads,
        ).filtered


        //
        // MIRNA BINDING SITES:
        //

        // Filtering miRNAs from ch_mature if they are not in ch_mirna_filtered.
        ch_uniq_mirnas = ch_mirna_filtered.map { _meta, path -> path }.splitCsv(sep: '\t').map { row -> row[0] }.unique().collect()

        ch_mature = ch_mature
            .map { _meta, path -> path }
            .splitFasta(record: [id: true, seqString: true])
            .combine(ch_uniq_mirnas.map { it -> [it] })
            .filter { record, _mirnas ->
                ch_uniq_mirnas.contains(record.id).value
            }
            .map { record, _mirnas ->
                ">${record.id}\n${record.seqString}"
            }
            .collectFile(name: 'mature_filtered.fa', newLine: true)
            .map { it -> [[id: 'mature_filtered'], it] }
    }

    MIRNA_BINDINGSITES(transcriptome_fasta, circrna_annotation, ch_mature, mirna_tools, mirna_min_tools)

    if (mirna_expression) {
        //
        // COMPUTE CORRELATION:
        //
        ch_binding_site_batches = MIRNA_BINDINGSITES.out.binding_sites
            .splitText(by: 100, file: true)
            .map { _meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file] }

        // Value channels, so that every binding site batch is combined with the same inputs
        COMPUTE_CORRELATIONS(ch_binding_site_batches, ch_mirna_filtered.first(), quantification_rds.first())

        ch_correlations = COMPUTE_CORRELATIONS.out.correlations
            .map { _meta, results -> results }
            .flatten()
            .collect()
            .map { results -> [[id: 'correlation'], results] }
    }

    emit:
    binding_sites = MIRNA_BINDINGSITES.out.binding_sites
    correlations  = ch_correlations
}
