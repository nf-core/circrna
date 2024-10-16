process CIRCRNA_FINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::circrna_finder=1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circrna_finder%3A1.2--pl5321hdfd78af_1' :
        'biocontainers/circrna_finder:1.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(star_input, stageAs: 'input/')

    output:
    tuple val(meta), path("${prefix}.filteredJunctions.bed"), emit: results
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 'v1.2'
    """
    postProcessStarAlignment.pl --starDir input/ --outDir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        cat: \$(cat --version | head -n 1 | sed -e 's/cat (GNU coreutils) //')
        circRNA_finder: $VERSION
    END_VERSIONS
    """
}
