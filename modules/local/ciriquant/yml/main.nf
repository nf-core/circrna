process CIRIQUANT_YML {
    label 'process_single'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'quay.io/biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    path gtf
    path fasta
    path bwa
    path hisat2

    output:
    path "travis.yml" , emit: yml

    when:
    task.ext.when == null || task.ext.when

    script:
    bwa_prefix = fasta.toString() == 'genome.fa' ? fasta.toString() : fasta.toString() - ~/.(fa|fasta)$/
    hisat2_prefix = fasta.toString() - ~/.(fa|fasta)$/
    fasta_path = fasta.toRealPath()
    gtf_path = gtf.toRealPath()
    bwa_path = bwa.toRealPath()
    hisat2_path = hisat2.toRealPath()
    """
    BWA=`which bwa`
    HISAT2=`which hisat2`
    STRINGTIE=`which stringtie`
    SAMTOOLS=`which samtools`

    touch travis.yml
    printf "name: ciriquant\ntools:\n  bwa: \$BWA\n  hisat2: \$HISAT2\n  stringtie: \$STRINGTIE\n  samtools: \$SAMTOOLS\n\nreference:\n  fasta: ${fasta_path}\n  gtf: ${gtf_path}\n  bwa_index: ${bwa_path}/${bwa_prefix}\n  hisat_index: ${hisat2_path}/${hisat2_prefix}" >> travis.yml
    """
}
