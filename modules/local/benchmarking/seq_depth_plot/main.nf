process SEQ_DEPTH_CORRELLATION {
    label "process_single"

    conda "bioconda::matplotlib=3.5.1 bioconda::seaborn=0.11.2"
    container 'uphl/seaborn'

    input:
        path(depth)
        tuple val(meta), path(bed)
    output:
       path("*.txt")  , optional:true, emit: report
    script:
        template "create_plots.py"
}