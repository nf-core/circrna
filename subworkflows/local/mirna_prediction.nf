// MODULES
include { DESEQ2_NORMALIZATION            } from '../../modules/local/deseq2/normalization'
include { MIRNA_FILTERING                 } from '../../modules/local/mirna_filtering'
include { COMPUTE_CORRELATIONS            } from '../../modules/local/compute_correlations'

// SUBWORKFLOWS
include { MIRNA_BINDINGSITES } from './mirna/mirna_bindingsites'

workflow MIRNA_PREDICTION {

    take:
    transcriptome_fasta
    circrna_bed12
    ch_mature
    ch_mirna
    transcript_counts
    quantification_rds

    main:
    ch_versions = Channel.empty()

    MIRNA_BINDINGSITES( transcriptome_fasta, circrna_bed12, ch_mature )
    ch_versions = ch_versions.mix(MIRNA_BINDINGSITES.out.versions)

    ADD_BACKSPLICE( circrna_fasta, [] )
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    //
    // MIRNA NORMALIZATION WORKFLOW:
    //

    ch_mirna_normalized = DESEQ2_NORMALIZATION( ch_mirna ).normalized

    ch_versions = ch_versions.mix(DESEQ2_NORMALIZATION.out.versions)

    ch_mirna_filtered = MIRNA_FILTERING( ch_mirna_normalized, 
                                         params.mirna_min_sample_percentage,  
                                         params.mirna_min_reads
                                         ).filtered
    
    ch_versions = ch_versions.mix(MIRNA_FILTERING.out.versions)

    //
    // COMPUTE CORREALTION:
    //

    //COMPUTE_CORRELATIONS( MIRNA_BINDINGSITES.out.binding_sites,
    //                      ch_mirna_filtered.map{meta, file -> file}.collect(),
    //                      transcript_counts.map{meta, file -> file}.collect()
    //                    )

    //ch_versions = ch_versions.mix(COMPUTE_CORRELATIONS.out.versions)
    
    emit:
    versions = ch_versions
}
