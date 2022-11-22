process TARGETSCAN_DATABASE {
    tag "$mature"
    label 'process_low'

    input:
    path(mature)

    output:
    path("mature.txt"), emit: mature_txt

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bash ${workflow.projectDir}/bin/targetscan_format.sh $mature
    """
}
