process CIRCRNA_FINDER_FILTER {
    tag "$meta.id"
    label 'process_low'

    container 'barryd237/circrna_finder'

    input:
    tuple val(meta), path(sam), path(junction), path(tab)
    path fasta
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_circrna_finder_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_circrna_finder.bed")      , emit: matrix
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

    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${prefix}.filteredJunctions.bed | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${prefix}_circrna_finder.bed

    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_circrna_finder.bed > ${prefix}_circrna_finder_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circRNA_finder: $VERSION
    END_VERSIONS
    """
}
