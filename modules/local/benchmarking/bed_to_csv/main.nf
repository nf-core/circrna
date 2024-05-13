process BED_TO_CSV {
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
        tuple val(id), path(bedfile1), path(bedfile2) from ch_joined
    output:
        path "${bedfile1.baseName}.csv"
        path "${bedfile2.baseName}.csv"
    script:

        """
        echo 'Process start'
        """
        template "process_bedfiles.py $bedfile1 ${bedfile1.baseName}.csv $bedfile2 ${bedfile2.baseName}.csv"
        

}