process OVERLAP_PLOT {
    label "process_single"

    conda "bioconda::seaborn=0.11.2"
    container 'uphl/seaborn'

    input:
        tuple val(meta), path(bed)
        path "versions.yml"

    output:
        path("*_mqc.png") , emit: plots
        path("versions.yml"), emit: versions

    script:
        template "create_plots.py"
}
