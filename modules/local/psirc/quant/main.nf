process PSIRC_QUANT {
    tag "${meta.id}"
    label 'process_high'

    container "registry.hub.docker.com/bigdatainbiomedicine/psirc"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${meta.id}/abundance.tsv"), emit: abundance_tsv
    tuple val(meta), path("${meta.id}/abundance.h5"), emit: abundance_h5
    tuple val(meta), path("${meta.id}"), emit: directory

    script:
    def single_end = meta.single_end ? "--single -l 76 -s 20" : ""
    """
    psirc-quant quant -t $task.cpus -i $index -o $meta.id $single_end $reads
    """
}
