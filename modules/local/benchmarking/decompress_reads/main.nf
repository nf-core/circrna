

process DECOMPRESS_READS {
    input:
    tuple val(meta), path(fastq_gz_1), path(fastq_gz_2)

    output:
    tuple val(meta), path("*.fq")

    script:
    """
    # Define output names by removing .gz and adding _decompressed.fq
    output_name1=\$(basename ${fastq_gz_1} .gz)_decompressed.fq
    output_name2=\$(basename ${fastq_gz_2} .gz)_decompressed.fq

    # Decompress the files
    gunzip -c ${fastq_gz_1} > \$output_name1
    gunzip -c ${fastq_gz_2} > \$output_name2
    """
}
