workflow STATISTICAL_TESTS {
    take:
    ch_quantification
    
    main:
    ch_versions = Channel.empty()

    ch_quantification.view()
    
    
}