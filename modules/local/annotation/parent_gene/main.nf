process PARENT_GENE {
    label 'process_high'

    conda "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0"
    container "registry.hub.docker.com/bigdatainbiomedicine/circ-annotation:latest"

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
        echo -e "\$chr\\t\$start\\t\$stop\\t\$sign" >> circs.bed
    done < IDs.txt

    mkdir -p bed12

    # Re-use annotation script to identify the host gene.
    parallel -j $task.cpus -a circs.bed annotate_outputs.sh $exon_boundary {}
    cat bed12/*.bed12.bed > master_bed12.bed.tmp
    awk 'BEGIN{FS=OFS="\\t"} {gsub(/,\$/,"",\$11);gsub(/,\$/,"",\$12)} 1' master_bed12.bed.tmp > master_bed12.bed && rm master_bed12.bed.tmp

    awk -v OFS="\\t" '{print \$4, \$14}' master_bed12.bed > circrna_host-gene.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        ucsc: $VERSION
    END_VERSIONS
    """
}
