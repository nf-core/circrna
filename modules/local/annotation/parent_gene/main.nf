process PARENT_GENE {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' :
        'quay.io/biocontainers/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' }"

    input:
    path circrna_matrix
    path gtf
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
    grep -vf ${workflow.projectDir}/bin/unwanted_biotypes.txt $gtf > filt.gtf

    # generate circrna BED file.
    tail -n +2 $circrna_matrix | awk '{print \$1}' > IDs.txt
    ID_to_BED.sh IDs.txt
    cat *.bed > merged.txt && rm IDs.txt && rm *.bed && mv merged.txt circs.bed

    # Re-use annotation script to identify the host gene.
    annotate_outputs.sh $exon_boundary &> annotation.log
    awk -v OFS="\t" '{print \$4, \$14}' master_bed12.bed > circrna_host-gene.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        ucsc: $VERSION
    END_VERSIONS
    """
}
