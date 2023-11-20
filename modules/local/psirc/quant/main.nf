process PSIRC_QUANT {
    tag "${meta.id}"
    label 'process_medium'

    container "registry.hub.docker.com/bigdatainbiomedicine/psirc"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("output/${meta.id}.tsv"), emit: abundance_tsv
    tuple val(meta), path("output/${meta.id}.h5"), emit: abundance_h5

    script:
    def single_end = meta.single_end ? "--single" : ""
    """
    psirc-quant quant -t $task.cpus -i $index -o output $single_end $reads

    mv output/abundance.tsv output/${meta.id}.tsv
    mv output/abundance.h5 output/${meta.id}.h5
    """
}