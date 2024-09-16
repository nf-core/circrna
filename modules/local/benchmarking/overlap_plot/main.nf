process OVERLAP_PLOT {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::seaborn=0.11.2"
    container 'community.wave.seqera.io/library/seaborn:0.13.2--ef0811a05c6fcc75'

    input:
        tuple val(meta), path(bed)

    output:
        path("*_mqc.png") , emit: plots
        path("versions.yml"), emit: versions

    script:
        template "create_plots.py"
}
