process CIRIQUANT {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/nicotru/ciriquant:35228e054dcad0b0'
        : 'docker.io/nicotru/ciriquant@sha256:5e59c7a31355f1814a1703e9cb27d347b7bb04871606a24c1acdacfc02c627f8'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(bwa)
    tuple val(meta6), path(hisat2)
    val strandedness

    output:
    tuple val(meta), path("${prefix}/${prefix}.gtf"), emit: gtf
    tuple val(meta), path("${prefix}/gene/${prefix}_genes.list"), emit: gene_list, optional: true
    tuple val(meta), path("${prefix}/gene/${prefix}_out.gtf"), emit: gene_gtf, optional: true
    tuple val("${task.process}"), val('bwa'), eval("bwa 2>&1 | sed -n 's/.*Version: //p' | head -n 1"), emit: versions_bwa, topic: versions
    tuple val("${task.process}"), val('ciriquant'), eval('CIRIquant --version 2>&1 | sed "s/CIRIquant //g"'), emit: versions_ciriquant, topic: versions
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version 2>&1 | sed -n '1{s/.* version //p}'"), emit: versions_hisat2, topic: versions
    tuple val("${task.process}"), val('pysam'), eval('python -c "import pysam; print(pysam.__version__)"'), emit: versions_pysam, topic: versions
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f2'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed -n '1{s/^samtools //p}'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version 2>&1'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def library_type = strandedness == 'auto' ? '' : strandedness == 'unstranded' ? '-l 0' : strandedness == 'forward' ? '-l 1' : '-l 2'
    def reads_string = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def bed_string = bed ? "--bed ${bed}" : ''
    """
    BWA=`which bwa`
    HISAT2=`which hisat2`
    STRINGTIE=`which stringtie`
    SAMTOOLS=`which samtools`

    BWA_FILE=`ls ${bwa}/*.bwt`
    BWA_PREFIX=`basename \$BWA_FILE .bwt`

    HISAT2_FILE=`ls ${hisat2}/*.1.ht2`
    HISAT2_PREFIX=`basename \$HISAT2_FILE .1.ht2`

    printf "name: ciriquant\\ntools:\\n  bwa: \$BWA\\n  hisat2: \$HISAT2\\n  stringtie: \$STRINGTIE\\n  samtools: \$SAMTOOLS\\n\\nreference:\\n  fasta: ${fasta}\\n  gtf: ${gtf}\\n  bwa_index: ${bwa}/\$BWA_PREFIX\\n  hisat_index: ${hisat2}/\$HISAT2_PREFIX" > config.yml

    CIRIquant \\
        -t ${task.cpus} \\
        ${reads_string} \\
        ${bed_string} \\
        --config config.yml \\
        -o ${prefix} \\
        -p ${prefix} \\
        ${library_type} \\
        ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/gene
    touch ${prefix}/${prefix}.gtf
    touch ${prefix}/gene/${prefix}_genes.list
    touch ${prefix}/gene/${prefix}_out.gtf

    """
}
