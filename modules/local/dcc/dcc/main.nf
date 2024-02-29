process DCC {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::circtools=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circtools:1.2.1--pyh7cba7a3_0' :
        'quay.io/biocontainers/circtools:1.2.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(pairs), path(mate1), path(mate2)
    path fasta
    path gtf

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if(meta.single_end){
        """
        sed -i 's/^chr//g' $gtf

        mkdir ${prefix} && mv ${prefix}.Chimeric.out.junction ${prefix} && printf "${prefix}/${prefix}.Chimeric.out.junction" > samplesheet

        DCC @samplesheet -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T ${task.cpus}

        awk '{print \$6}' CircCoordinates >> strand
        paste CircRNACount strand | tail -n +2 | awk -v OFS="\\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${prefix}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dcc: \$(DCC --version)
        END_VERSIONS
        """
    }else{
        """
        sed -i 's/^chr//g' $gtf

        mkdir ${prefix} && mv ${prefix}.Chimeric.out.junction ${prefix} && printf "${prefix}/${prefix}.Chimeric.out.junction" > samplesheet
        mkdir ${prefix}_mate1 && mv ${prefix}_mate1.Chimeric.out.junction ${prefix}_mate1 && printf "${prefix}_mate1/${prefix}_mate1.Chimeric.out.junction" > mate1file
        mkdir ${prefix}_mate2 && mv ${prefix}_mate2.Chimeric.out.junction ${prefix}_mate2 && printf "${prefix}_mate2/${prefix}_mate2.Chimeric.out.junction" > mate2file

        DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T ${task.cpus}

        awk '{print \$6}' CircCoordinates >> strand
        paste CircRNACount strand | tail -n +2 | awk -v OFS="\\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${prefix}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dcc: \$(DCC --version)
        END_VERSIONS
        """
    }
}
