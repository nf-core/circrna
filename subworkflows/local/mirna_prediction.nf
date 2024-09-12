// MODULES
include { BIOAWK as ADD_BACKSPLICE                } from '../../modules/nf-core/bioawk'
include { MIRNA_NORMALIZATION                     } from '../../modules/local/deseq2/mirna_normalization'
include { TX_NORMALIZATION                        } from '../../modules/local/deseq2/tx_normalization'
include { MIRNA_FILTERING                         } from '../../modules/local/mirna_filtering'
include { COMPUTE_CORRELATIONS                    } from '../../modules/local/compute_correlations'
include { SPONGE                                  } from '../../modules/local/sponge'
include { SPONGE_EFFECTS                          } from '../../modules/local/sponge_effects'

// SUBWORKFLOWS
include { MIRNA_BINDINGSITES } from './mirna/mirna_bindingsites'

workflow MIRNA_PREDICTION {

    take:
    transcriptome_fasta
    circrna_annotation
    ch_mature
    ch_mirna
    quantification_rds
    tx_counts

    main:
    ch_versions = Channel.empty()

    //
    // MIRNA NORMALIZATION WORKFLOW:
    //

    if (params.mirna_expression) {

        ch_mirna_normalized = MIRNA_NORMALIZATION(ch_mirna).normalized

        ch_versions = ch_versions.mix(MIRNA_NORMALIZATION.out.versions)

        ch_mirna_filtered = MIRNA_FILTERING(ch_mirna_normalized,
                                            params.mirna_min_sample_percentage,
                                            params.mirna_min_reads
                                            ).filtered

        ch_versions = ch_versions.mix(MIRNA_FILTERING.out.versions)

        //
        // MIRNA BINDING SITES:
        //

        // Filtering miRNAs from ch_mature if they are not in ch_mirna_filtered.
        ch_uniq_mirnas = ch_mirna_filtered.map{ meta, path -> path }.splitCsv( sep: '\t' ).map{ it[0] }.unique().collect()

        ch_mature = ch_mature
            .map{ meta, path ->
                path
            }
            .splitFasta( record: [id:true, seqString:true] )
            .combine(ch_uniq_mirnas.map{ it -> [it]}) // Not sure why this mapping is necessary but I think it is
            .filter{ record, mirnas ->
                ch_uniq_mirnas.contains(record.id).value
            }.map{ record, mirnas ->
                ">${record.id}\n${record.seqString}"
            }
            .collectFile( name: 'mature_filtered.fa', newLine: true)
            .map{ it -> [[id: 'mature_filtered'], it]}
    }

    MIRNA_BINDINGSITES(transcriptome_fasta, circrna_annotation, ch_mature)
    ch_versions = ch_versions.mix(MIRNA_BINDINGSITES.out.versions)

    if (params.mirna_expression) {
        //
        // COMPUTE CORRELATION:
        //
        ch_binding_site_batches = MIRNA_BINDINGSITES.out.targets
            .splitText(by: 100, file: true)
            .map{ meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file]}

        COMPUTE_CORRELATIONS(ch_binding_site_batches, ch_mirna_filtered, quantification_rds)

        ch_correlation_results = COMPUTE_CORRELATIONS.out.correlations
            .map{meta, results -> results}
            .flatten().collect()
            .map{results -> [[id: 'correlation'], results]}

        ch_versions = ch_versions.mix(COMPUTE_CORRELATIONS.out.versions)
    }

    // if (params.sponge)


    ch_tx_normalized = TX_NORMALIZATION(tx_counts).normalized

    ch_versions = ch_versions.mix(TX_NORMALIZATION.out.versions)
    
    SPONGE(MIRNA_BINDINGSITES.out.binding_sites, ch_tx_normalized, ch_mirna_filtered)

    ch_versions = ch_versions.mix(SPONGE.out.versions)

    SPONGE_EFFECTS(SPONGE.out.sponge_data)
    

    emit:
    versions = ch_versions
}
