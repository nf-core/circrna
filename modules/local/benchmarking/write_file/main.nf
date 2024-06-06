process WRITE {
    label "process_single"

    conda "bioconda::matplotlib=3.5.1 bioconda::seaborn=0.11.2"
    container 'uphl/seaborn'

    input:
        tuple val(meta), path(value)
    output:
       path("*.txt")  , optional:true, emit: report
    script:
        template "write.py"
}