process CIRCTOOLS_DETECT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f523ec9d03407da839bae40c76b9e49725d8b5273e9e167f7b5e081392f690d/data'
        : 'community.wave.seqera.io/library/circtools:2.0.5.post2--20ffd704d69916e2'}"

    input:
    tuple val(meta), path(paired, stageAs: 'paired.junctions'), path(mate1, stageAs: 'mate1.junctions'), path(mate2, stageAs: 'mate2.junctions')
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    val strandedness

    output:
    tuple val(meta), path("${prefix}_reads.junctions"), emit: reads
    tuple val(meta), path("${prefix}_coordinates.tsv"), emit: coordinates
    tuple val(meta), path("${prefix}_counts.tsv"), emit: counts

    tuple val("${task.process}"), val('circtools'), eval('circtools -V'), emit: versions_circtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def strand_args = strandedness == 'auto' || strandedness == 'unstranded' ? '-N' : strandedness == 'forward' ? '' : '-ss'

    def matefile_commands = meta.single_end
        ? ''
        : "printf ${mate1} > mate1file && printf ${mate2} > mate2file"

    def mate_args = meta.single_end ? '' : '-mt1 @mate1file -mt2 @mate2file -Pi'

    // circtools only writes tmp_printcirclines in paired-end independent mode (-Pi). For single-end data the
    // equivalent is built from the STAR chimeric junctions: no header or comment lines, 14 STAR columns plus the sample.
    def reads_command = meta.single_end
        ? """awk -F'\\t' -v OFS='\\t' -v sample='${meta.id}' '!/^#/ && \$1 != "chr_donorA" && NF >= 14 { for (i = 1; i <= 14; i++) printf "%s\\t", \$i; print sample }' paired.junctions > ${prefix}_reads.junctions"""
        : "mv _tmp_circtools/tmp_printcirclines.[0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z] ${prefix}_reads.junctions"
    """
    printf "${paired}" > samplesheet
    ${matefile_commands}

    circtools detect @samplesheet ${mate_args} -D -an ${gtf} ${args} -F -M -k -Nr 1 1 -A ${fasta} ${strand_args} -T ${task.cpus}

    ${reads_command}
    mv CircCoordinates ${prefix}_coordinates.tsv
    mv CircRNACount ${prefix}_counts.tsv

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_reads.junctions
    touch ${prefix}_coordinates.tsv
    touch ${prefix}_counts.tsv

    """
}
