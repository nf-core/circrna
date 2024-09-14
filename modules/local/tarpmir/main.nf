process TARPMIR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/miranda_viennarna_pip_numpy_pruned:d46fd1a59195bfaf' :
        'community.wave.seqera.io/library/miranda_viennarna_pip_numpy_pruned:f4624328aba4ca2b' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(mature)

    output:
    tuple val(meta), path("${prefix}.txt")           , emit: bindings
    tuple val(meta), path("*.txt")                   , emit: txt
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'TarPmiR_threading.py'
    
    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miranda: \$(echo \$(miranda -v | sed -n 4p | sed 's/^.*miranda v//; s/microRNA.*\$//' ))
        RNAplfold: \$(echo \$(RNAplfold -V | sed 's/^.*RNAplfold //' ))
        RNAduplex: \$(echo \$(RNAduplex -V | sed 's/^.*RNAduplex //' ))
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
        pymysql: \$(python -c "import pymysql; print(pymysql.__version__)")
    END_VERSIONS
    """
}
