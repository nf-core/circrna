process FAIL_ON_EMPTY {
    tag "$meta.id"

    input:
    tuple val(meta), path(bed)
    path(waitFor, stageAs: 'waitFor/*.txt')
    val(message)

    exec:
    if (!bed) {
        log.error message

        exit 1
    }

    stub:
    if (!bed) {
        log.error message

        exit 1
    }
}
