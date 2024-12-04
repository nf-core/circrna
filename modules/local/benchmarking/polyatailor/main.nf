process POLYATAILOR {
    tag "$meta.id"
    label "process_high"

    container 'docker://faboehm/polyatailor-env'

    input:
        tuple val(meta),
        path(fastq)
        path(bam)

    output:
        path("*.tsv"), emit: tails
        path("versions.yml"), emit: versions

    script:
    template'polyatailor.R'
    
}

