process PSIRC_QUANT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/psirc:1.0.0--he1fd2f9_0' :
        'biocontainers/psirc:1.0.0--he1fd2f9_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(chrom_sizes)
    val(bootstrap_samples)

    output:
    tuple val(meta), path("${meta.id}"), emit: directory
    path "versions.yml"                , emit: versions

    script:
    def single_end = meta.single_end ? "--single -l 76 -s 20" : ""
    def genomebam = gtf ? "--genomebam -g $gtf" : ""
    def chromosomes = chrom_sizes ? "-c $chrom_sizes" : ""
    """
    psirc-quant quant -t $task.cpus -i $index -o $meta.id $single_end $reads -b $bootstrap_samples $genomebam $chromosomes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc-quant: \$(psirc-quant version | sed -n 's/^psirc-quant, version \\([0-9.]*\\).*\$/\\1/p')
    END_VERSIONS
    """
}
