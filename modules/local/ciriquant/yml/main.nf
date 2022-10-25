process CIRIQUANT_YML {
    label 'process_single'

    container 'barryd237/ciriquant_v1.0.1'

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
    BWA=`whereis bwa | cut -f2 -d':'`
    HISAT2=`whereis hisat2 | cut -f2 -d':'`
    STRINGTIE=`whereis stringtie | cut -f2 -d':'`
    SAMTOOLS=`whereis samtools | cut -f2 -d':' | awk '{print \$1}'`

    touch travis.yml
    printf "name: ciriquant\ntools:\n  bwa:\$BWA\n  hisat2:\$HISAT2\n  stringtie:\$STRINGTIE\n  samtools: \$SAMTOOLS\n\nreference:\n  fasta: ${fasta_path}\n  gtf: ${gtf_path}\n  bwa_index: ${bwa_path}/${bwa_prefix}\n  hisat_index: ${hisat2_path}/${hisat2_prefix}" >> travis.yml
    """
}
