process PSIRC_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/psirc:1.0.0--he1fd2f9_0' :
        'biocontainers/psirc:1.0.0--he1fd2f9_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("psirc.index"), emit: index
    path "versions.yml",                  emit: versions

    script:
    """
    psirc-quant index -i psirc.index --make-unique $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc-quant: \$(psirc-quant version | sed -n 's/^psirc-quant, version \\([0-9.]*\\).*\$/\\1/p')
    END_VERSIONS
    """
}
