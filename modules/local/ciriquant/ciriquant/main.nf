process CIRIQUANT {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nicotru/ciriquant:1.0.1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(bwa)
    tuple val(meta6), path(hisat2)

    output:
    tuple val(meta), path("${prefix}/${prefix}.gtf")            , emit: gtf
    tuple val(meta), path("${prefix}/gene/${prefix}_genes.list"), emit: gene_list, optional: true
    tuple val(meta), path("${prefix}/gene/${prefix}_out.gtf") , emit: gene_gtf, optional: true
    path("${prefix}")                                           , emit: results
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.0'
    def strandedness = meta.strandedness ?: 'auto'
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        stringtie: \$(stringtie --version 2>&1)
        hisat2: $VERSION
    END_VERSIONS
    """
}
