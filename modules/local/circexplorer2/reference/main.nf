process CIRCEXPLORER2_REFERENCE {
    tag "$gtf"
    label 'process_single'

    conda "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' :
        'quay.io/biocontainers/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' }"

    input:
    path gtf

    output:
    path("${prefix}.txt"), emit: txt
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = gtf.simpleName
    def VERSION = '377'
    """
    gtfToGenePred \
        $args \
        $gtf \
        ${prefix}.genepred

    awk -v OFS="\\t" '{print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' ${prefix}.genepred > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        ucsc: $VERSION
    END_VERSIONS
    """
}
