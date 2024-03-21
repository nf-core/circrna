process CIRIQUANT_YML {
    label 'process_single'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'quay.io/biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(fasta)
    path bwa
    path hisat2

    output:
    path "travis.yml" , emit: yml

    when:
    task.ext.when == null || task.ext.when

    script:
    hisat2_prefix = meta2.id
    fasta_path = fasta.toRealPath()
    gtf_path = gtf.toRealPath()
    bwa_path = bwa.toRealPath()
    hisat2_path = hisat2.toRealPath()
    """
    BWA=`which bwa`
    HISAT2=`which hisat2`
    STRINGTIE=`which stringtie`
    SAMTOOLS=`which samtools`

    # Get first file in bwa index directory
    BWA_FILE=`ls ${bwa_path}/*.bwt`
    BWA_FILE=`basename \$BWA_FILE .bwt`

    touch travis.yml
    printf "name: ciriquant\ntools:\n  bwa: \$BWA\n  hisat2: \$HISAT2\n  stringtie: \$STRINGTIE\n  samtools: \$SAMTOOLS\n\nreference:\n  fasta: ${fasta_path}\n  gtf: ${gtf_path}\n  bwa_index: ${bwa_path}/\$BWA_FILE\n  hisat_index: ${hisat2_path}/${hisat2_prefix}" >> travis.yml
    """
}
