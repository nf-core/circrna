process STRINGTIE_PREPDE {
    label 'process_low'

    conda "bioconda::stringtie=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    path gtf

    output:
    path "transcript_count_matrix.csv" , emit: transcript_matrix
    path "gene_count_matrix.csv"       , emit: gene_matrix
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for file in \$(ls *.gtf); do sample_id=\${file%".transcripts.gtf"}; touch samples.txt; printf "\$sample_id\\t\$file\\n" >> samples.txt ; done

    prepDE.py -i samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
