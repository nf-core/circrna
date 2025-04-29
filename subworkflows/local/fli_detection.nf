workflow FLI_DETECTION {
    take:
    ch_bsj_bed
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    

    emit:
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
