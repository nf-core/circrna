process RRNA_CORRELATION {
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
        val bed_real_list   // Collected list of tuples for real BEDs
        val bed_bench_list  // Collected list of tuples for benchmarking BEDs
        val rRNA_real_list  // Collected list of tuples for real rRNA summaries
        val rRNA_bench_list // Collected list of tuples for benchmarking rRNA summaries

    script:
        template "rRNA_corr.py"
}
