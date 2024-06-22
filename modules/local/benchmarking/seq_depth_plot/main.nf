process SEQ_DEPTH_CORRELLATION {
    label "process_single"

    conda "bioconda::matplotlib=3.5.1 bioconda::seaborn=0.11.2"
    container 'https://depot.galaxyproject.org/singularity/bionumpy:0.2.12--pyha8f3691_0'

    input:
        tuple val(meta),
        path(bed)
        tuple val(depth_meta),
        path(depth)
    output:
       path("*.txt")  , emit: report
    script:
        template "create_plots.py"
}