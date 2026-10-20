process STRINGTIE_STRINGTIE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56eda015168efa0a6b6c814ae70827d2a05906d8f6e21450d6ccad05b29c7437/data'
        : 'community.wave.seqera.io/library/stringtie:3.0.3--f9bc4844393ae0fa'}"

    input:
    tuple val(meta), path(bam)
    path  annotation_gtf

    output:
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), path("*.abundance.txt")  , emit: abundance
    tuple val(meta), path("*.coverage.gtf")   , optional: true, emit: coverage_gtf
    tuple val(meta), path("*.ballgown/*.ctab"), optional: true, emit: ballgown
    tuple val("${task.process}"), val('stringtie'), eval("stringtie --version"), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    def ballgown  = annotation_gtf ? "-b ${prefix}.ballgown" : ""
    def coverage  = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandedness \\
        $reference \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene.abundance.txt \\
        $coverage \\
        $ballgown \\
        -p $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.gene.abundance.txt
    touch ${prefix}.coverage.gtf
    mkdir -p ${prefix}.ballgown
    touch ${prefix}.ballgown/e_data.ctab ${prefix}.ballgown/i_data.ctab ${prefix}.ballgown/t_data.ctab ${prefix}.ballgown/e2t.ctab ${prefix}.ballgown/i2t.ctab
    """
}
