process PLOT_POLYATAILS {
    label "process_single"

    conda "bioconda::seaborn=0.11.2"
    container 'community.wave.seqera.io/library/seaborn:0.13.2--ef0811a05c6fcc75'

    input:
        path(real)
        path(benchmarking)
    output:
        path("*.tsv")  , emit: average_tails
    script:
        template "create_plots.py"
}