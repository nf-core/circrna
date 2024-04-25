process ANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(intersection)
    val(exon_boundary)

    output:
    tuple val(meta), path("${out_bed}"), emit: bed
    tuple val(meta), path("${out_gtf}"), emit: gtf
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}.annotation"
    out_bed = "${prefix}.bed"
    out_gtf = "${prefix}.gtf"
    """
    annotation.py --input ${intersection} --exon_boundary ${exon_boundary} --output_bed ${out_bed} --output_gtf ${out_gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
