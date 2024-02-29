process CIRCRNA_FINDER_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::circrna_finder=1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circrna_finder%3A1.2--pl5321hdfd78af_1' :
        'quay.io/biocontainers/circrna_finder:1.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(sam), path(junction), path(tab)
    path fasta
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_circrna_finder_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_circrna_finder.bed")      , emit: matrix
    tuple val(meta), path("*filteredJunctions*")               , emit: intermediates
    path  "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 'v1.2'
    """
    mkdir -p star_dir && mv *.tab *.junction *.sam star_dir
    postProcessStarAlignment.pl --starDir star_dir/ --outDir ./

    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${prefix}.filteredJunctions.bed | awk  -v OFS="\\t" -F"\\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${prefix}_circrna_finder.bed

    awk -v OFS="\\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_circrna_finder.bed > ${prefix}_circrna_finder_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        cat: \$(cat --version | head -n 1 | sed -e 's/cat (GNU coreutils) //')
        circRNA_finder: $VERSION
    END_VERSIONS
    """
}
