process FASTA {
    tag "${meta.id}:${meta.tool}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("${prefix}.fa"), emit: analysis_fasta
    path("${prefix}.fasta")              , emit: publish_fasta
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377'
    """
    ## FASTA sequences (bedtools does not like the extra annotation info - split will not work properly)
    cut -d\$'\t' -f1-12 ${base}.bed > bed12.tmp
    bedtools getfasta -fi $fasta -bed bed12.tmp -s -split -name > circ_seq.tmp

    ## clean fasta header
    grep -A 1 '>' circ_seq.tmp | cut -d: -f1,2,3 > ${base}.fa && rm circ_seq.tmp

    ## add backsplice sequence for miRanda Targetscan, publish canonical FASTA to results.
    rm $fasta
    bash ${workflow.projectDir}/bin/backsplice_gen.sh ${base}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}


