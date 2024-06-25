process LOCATION_PLOT {
    label "process_single"

    conda "bioconda::matplotlib=3.5.1 bioconda::seaborn=0.11.2"
    container 'uphl/seaborn'

    input:
        tuple val(id), path(bedfile1), path(bedfile2)
    output:
        path("*_mqc.png") , emit: report
    script:
        template "create_plots.py"
}
