process PARENT_GENE {
    label 'process_high'

    conda "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' :
        'quay.io/biocontainers/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' }"

    input:
    path circrna_matrix
    path gtf
    path biotypes
    val exon_boundary

    output:
    path "circrna_host-gene.txt" , emit: circ_host_map
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '377'
    """
    # remove redundant biotypes from GTF.
    grep -vf $biotypes $gtf > filt.gtf

    # generate circrna BED file.
    tail -n +2 $circrna_matrix | awk '{print \$1}' > IDs.txt

    # ID_to_BED
    while IFS='' read -r line; do
        name=\$(echo \$line)
        chr=\$(echo \$line | cut -d: -f1)
        start=\$(echo \$line | cut -d- -f1 | cut -d: -f2)
        stop=\$(echo \$line | cut -d- -f2 | cut -d: -f1)
        sign=\$(echo \$line | cut -d: -f3)
        echo -e "\$chr\\t\$start\\t\$stop\\t\$name\\t0\\t\$sign" >> \${name}.bed
    done < IDs.txt

    cat *.bed > merged.txt && rm IDs.txt && rm *.bed && mv merged.txt circs.bed

    # Re-use annotation script to identify the host gene.
    annotate_outputs.sh $exon_boundary &> annotation.log
    awk -v OFS="\\t" '{print \$4, \$14}' master_bed12.bed > circrna_host-gene.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        ucsc: $VERSION
    END_VERSIONS
    """
}
