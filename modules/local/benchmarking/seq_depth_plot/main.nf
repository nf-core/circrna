process SEQ_DEPTH_CORRELLATION {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::matplotlib=3.5.1 bioconda::seaborn=0.11.2"
    container 'https://depot.galaxyproject.org/singularity/bionumpy:0.2.12--pyha8f3691_0'

    input:
        tuple val(meta),
        path(bed)
        path(depth)
    output:
        path("*.tsv")  , emit: report
        path("versions.yml"), emit: versions
    script:
        template "correlation.py"
}
