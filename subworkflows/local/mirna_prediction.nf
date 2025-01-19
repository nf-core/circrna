// MODULES
include { BIOAWK as ADD_BACKSPLICE                } from '../../modules/nf-core/bioawk'
include { MIRNA_NORMALIZATION                     } from '../../modules/local/deseq2/mirna_normalization'
include { GENE_NORMALIZATION                      } from '../../modules/local/deseq2/gene_normalization'
include { MIRNA_FILTERING                         } from '../../modules/local/mirna/filtering'
include { MIRNA_COMPUTECORRELATIONS               } from '../../modules/local/mirna/computecorrelations'
include { SPONGE_SPONGE                           } from '../../modules/local/sponge/sponge'
include { SPONGE_SPONGEFFECTS                     } from '../../modules/local/sponge/spongeffects'

// SUBWORKFLOWS
include { MIRNA_BINDINGSITES } from './mirna/mirna_bindingsites'

workflow MIRNA_PREDICTION {

    take:
    transcriptome_fasta
    circrna_annotation
    ch_mature
    ch_mirna
    tx_counts
    quantification_rds

    main:
    ch_versions = Channel.empty()

    //
    // MIRNA NORMALIZATION WORKFLOW
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
        // MIRNA BINDING SITES
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
        // COMPUTE CORRELATION
        //
        ch_binding_site_batches = MIRNA_BINDINGSITES.out.targets
            .splitText(by: 100, file: true)
            .map{ meta, file -> [[id: "batch_" + file.baseName.split("\\.").last()], file]}

        MIRNA_COMPUTECORRELATIONS(ch_binding_site_batches, ch_mirna_filtered, quantification_rds)

        ch_correlation_results = MIRNA_COMPUTECORRELATIONS.out.correlations
            .map{meta, results -> results}
            .flatten().collect()
            .map{results -> [[id: 'correlation'], results]}

        ch_versions = ch_versions.mix(MIRNA_COMPUTECORRELATIONS.out.versions)

        //
        // SPONGE
        //
        ch_gene_normalized = GENE_NORMALIZATION(tx_counts).normalized
        ch_versions = ch_versions.mix(GENE_NORMALIZATION.out.versions)

        SPONGE(MIRNA_BINDINGSITES.out.binding_sites, ch_gene_normalized, ch_mirna_filtered)
        ch_versions = ch_versions.mix(SPONGE.out.versions)

        SPONGE_EFFECTS(SPONGE.out.sponge_data)
        ch_versions = ch_versions.mix(SPONGE_EFFECTS.out.versions)
    }

    emit:
    versions = ch_versions
}
