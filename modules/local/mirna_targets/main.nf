process MIRNA_TARGETS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/36dcc9054fa4a4d366be5665430b8d3c61b9eb949653fafbc243d055f627bba4/data'
        : 'community.wave.seqera.io/library/mirna_targets:f0ede764f11e5da4'}"

    input:
    tuple val(meta), path(targetscan), path(miranda)

    output:
    tuple val(meta), path("${prefix}.mirna_targets.txt"), emit: results
    tuple val("${task.process}"), val('awk'), eval('awk --version | head -n1 | cut -d\' \' -f3 | sed \'s/,//g\''), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed "s/bedtools v//"'), emit: versions_bedtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## reformat and sort miRanda, TargetScan outputs, convert to BED for overlaps.
    tail -n +2 $targetscan | sort -k1,1 -k4n | awk -v OFS="\\t" '{print \$1, \$2, \$4, \$5, \$9}' | awk -v OFS="\\t" '{print \$2, \$3, \$4, \$1, "0", \$5}' > targetscan.bed
    tail -n +2 $miranda | sort -k2,2 -k7n | awk -v OFS="\\t" '{print \$2, \$1, \$3, \$4, \$7, \$8}' | awk -v OFS="\\t" '{print \$2, \$5, \$6, \$1, \$3, \$4}' > miranda.bed

    ## intersect, consolidate miRanda, TargetScan information about miRs.
    ## -wa to output miRanda hits - targetscan makes it difficult to resolve duplicate miRNAs at MRE sites.
    bedtools intersect -a miranda.bed -b targetscan.bed -wa > ${prefix}.mirnas.tmp
    bedtools intersect -a targetscan.bed -b miranda.bed | awk '{print \$6}' > mirna_type

    ## remove duplicate miRNA entries at MRE sites.
    ## strategy: sory by circs, sort by start position, sort by site type - the goal is to take the best site type (i.e rank site type found at MRE site).
    paste ${prefix}.mirnas.tmp mirna_type | sort -k3n -k2n -k7r | awk -v OFS="\\t" '{print \$4,\$1,\$2,\$3,\$5,\$6,\$7}' | awk -F "\\t" '{if (!seen[\$1,\$2,\$3,\$4,\$5,\$6]++)print}' | sort -k1,1 -k3n > ${prefix}.mirna_targets.tmp
    echo -e "circRNA\\tmiRNA\\tStart\\tEnd\\tScore\\tEnergy_KcalMol\\tSite_type" | cat - ${prefix}.mirna_targets.tmp > ${prefix}.mirna_targets.txt

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mirna_targets.txt

    """
}
