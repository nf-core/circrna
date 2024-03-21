process MIRNA_TARGETS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(targetscan), path(miranda), path(bed12)

    output:
    tuple val(meta), path("${prefix}.mirna_targets.txt"), emit: results
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## reformat and sort miRanda, TargetScan outputs, convert to BED for overlaps.
    tail -n +2 $targetscan | sort -k1,1 -k4n | awk -v OFS="\\t" '{print \$1, \$2, \$4, \$5, \$9}' | awk -v OFS="\t" '{print \$2, \$3, \$4, \$1, "0", \$5}' > targetscan.bed
    tail -n +2 $miranda | sort -k2,2 -k7n | awk -v OFS="\\t" '{print \$2, \$1, \$3, \$4, \$7, \$8}' | awk -v OFS="\t" '{print \$2, \$5, \$6, \$1, \$3, \$4}' | sed 's/^[^-]*-//g' > miranda.bed

    ## intersect, consolidate miRanda, TargetScan information about miRs.
    ## -wa to output miRanda hits - targetscan makes it difficult to resolve duplicate miRNAs at MRE sites.
    bedtools intersect -a miranda.bed -b targetscan.bed -wa > ${prefix}.mirnas.tmp
    bedtools intersect -a targetscan.bed -b miranda.bed | awk '{print \$6}' > mirna_type

    ## remove duplicate miRNA entries at MRE sites.
    ## strategy: sory by circs, sort by start position, sort by site type - the goal is to take the best site type (i.e rank site type found at MRE site).
    paste ${prefix}.mirnas.tmp mirna_type | sort -k3,3 -k2n -k7r | awk -v OFS="\\t" '{print \$4,\$1,\$2,\$3,\$5,\$6,\$7}' | awk -F "\\t" '{if (!seen[\$1,\$2,\$3,\$4,\$5,\$6]++)print}' | sort -k1,1 -k3n > ${prefix}.mirna_targets.tmp
    echo -e "circRNA\\tmiRNA\\tStart\\tEnd\\tScore\\tEnergy_KcalMol\\tSite_type" | cat - ${prefix}.mirna_targets.tmp > ${prefix}.mirna_targets.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
