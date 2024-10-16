process FAIL_ON_EMPTY {
    tag "$meta.id"

    input:
    tuple val(meta), path(bed)
    path(waitFor, stageAs: 'waitFor*.txt')

    exec:
    if (!bed) {
        log.error ((
            "No circular RNAs were found by at least ${params.min_tools} tools and in at least ${params.min_samples} samples.\n") +
            "Feel free to check the preliminary results in '${params.outdir}'\n" +
            (params.save_intermediates ? "" :
            "You can enable saving intermediate files by setting the parameter 'save_intermediates' to 'true'."))

        exit 1
    }
}
