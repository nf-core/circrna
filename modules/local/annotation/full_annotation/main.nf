process ANNOTATION {
    tag "${meta.id}:${meta.tool}"
    label 'process_high'

    conda "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0 conda-forge::parallel=20230922"
    container "registry.hub.docker.com/bigdatainbiomedicine/circ-annotation:latest"

    input:
    tuple val(meta), path(bed)
    path gtf
    path biotypes
    val exon_boundary


    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377'
    """
    grep -vf $biotypes $gtf > filt.gtf
    mv $bed circs.bed

    mkdir -p bed12

    parallel -j $task.cpus -a circs.bed annotate_outputs.sh $exon_boundary {}
    cat bed12/*.bed12.bed > master_bed12.bed.tmp
    awk 'BEGIN{FS=OFS="\\t"} {gsub(/,\$/,"",\$11);gsub(/,\$/,"",\$12)} 1' master_bed12.bed.tmp > master_bed12.bed && rm master_bed12.bed.tmp

    mv master_bed12.bed ${prefix}.bed.tmp

    awk -v FS="\\t" '{print \$11}' ${prefix}.bed.tmp > mature_len.tmp
    awk -v FS="," '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' mature_len.tmp > mature_length

    paste ${prefix}.bed.tmp mature_length > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        ucsc: $VERSION
    END_VERSIONS
    """
}
