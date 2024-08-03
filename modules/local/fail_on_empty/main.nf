process FAIL_ON_EMPTY {
    tag "$meta.id"

    input:
    tuple val(meta), path(bed)
    path(waitFor, stageAs: 'waitFor*.txt')

    exec:
    if (!bed) {
        if (params.tool_filter <= 1) {
            log.error "No circular RNAs were found by any tool in any sample.\n" +
                "Please check the input data."
        } else {
            log.error "No circular RNAs were found by at least ${params.tool_filter} tools in any sample.\n" +
                "Please check the input data or lower the tool_filter parameter."
        }
        exit 1
    }
}
