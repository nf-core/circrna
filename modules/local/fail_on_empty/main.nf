process FAIL_ON_EMPTY {
    tag "$meta.id"

    input:
    tuple val(meta), path(bed)
    path(waitFor, stageAs: 'waitFor*.txt')

    exec:
    if (!bed) {
        log.error ((params.tool_filter <= 1 ?
            "No circular RNAs were found by any tool in any sample.\n" :
            "No circular RNAs were found by at least ${params.tool_filter} tools in any sample.\n") +
            "Feel free to check the preliminary results in '${params.outdir}'\n" +
            (params.save_intermediates ? "" :
            "You can enable saving intermediate files by setting the parameter 'save_intermediates' to 'true'."))

        exit 1
    }
}
