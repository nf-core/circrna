include { CIRI3_CIRI3 as CIRI3 } from '../../../../modules/local/ciri3/ciri3'
include { JOIN_SAMPLES as JOIN_BSJ } from '../../../../modules/local/matrix/join_samples'
include { JOIN_SAMPLES as JOIN_FSJ } from '../../../../modules/local/matrix/join_samples'

workflow CIRI3_QUANT {
    take:
    sam
    bed12
    ch_fasta
    ch_gtf

    main:
    ch_versions = channel.empty()

    // CIRI3 -UC takes chromosome, 1-based start, end. BED12 starts are 0-based.
    CIRI3_UC_LIST(bed12)
    ch_versions = ch_versions.mix(CIRI3_UC_LIST.out.versions)

    CIRI3(sam, ch_fasta, ch_gtf, CIRI3_UC_LIST.out.list.map { meta, list -> list })
    ch_versions = ch_versions.mix(CIRI3.out.versions)

    CIRI3_UC_MATRIX(CIRI3.out.txt, CIRI3_UC_LIST.out.list)
    ch_versions = ch_versions.mix(CIRI3_UC_MATRIX.out.versions)

    JOIN_BSJ(
        CIRI3_UC_MATRIX.out.bsj.map { meta, table -> [[id: 'ciri3_bsj'], meta.id, table] }.groupTuple(),
        'circRNA_ID'
    )
    ch_versions = ch_versions.mix(JOIN_BSJ.out.versions)

    JOIN_FSJ(
        CIRI3_UC_MATRIX.out.fsj.map { meta, table -> [[id: 'ciri3_fsj'], meta.id, table] }.groupTuple(),
        'circRNA_ID'
    )
    ch_versions = ch_versions.mix(JOIN_FSJ.out.versions)

    emit:
    bsj      = JOIN_BSJ.out.joined.map { meta, table -> [[id: 'ciri3_bsj', quantification: 'ciri3'], table] }
    fsj      = JOIN_FSJ.out.joined.map { meta, table -> [[id: 'ciri3_fsj', quantification: 'ciri3'], table] }
    versions = ch_versions
}

process CIRI3_UC_LIST {
    tag "consensus"
    label 'process_single'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path('consensus.ciri3.list'), emit: list
    path 'versions.yml', emit: versions

    script:
    """
    awk 'BEGIN { OFS="\\t" } !/^#/ && NF >= 3 { if (\$2 < 0) exit 2; print \$1, \$2 + 1, \$3 }' ${bed} > consensus.ciri3.list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>&1 | sed -n '1s/.*Awk //; s/,.*//p')
    END_VERSIONS
    """

    stub:
    """
    touch consensus.ciri3.list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}

process CIRI3_UC_MATRIX {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(result)
    tuple val(meta2), path(circ_list)

    output:
    tuple val(meta), path("${meta.id}.bsj.tsv"), emit: bsj
    tuple val(meta), path("${meta.id}.fsj.tsv"), emit: fsj
    path 'versions.yml', emit: versions

    script:
    """
    # UC matrices contain only observed rows. Re-key by exact genomic coordinates
    # and iterate the supplied list to retain the complete consensus universe.
    awk 'BEGIN { FS=OFS="\\t" } NR==FNR { key=\$1 ":" \$2 "|" \$3; ids[++n]=key; chr[key]=\$1; start[key]=\$2; end[key]=\$3; next } FNR>1 && NF>=4 { key=\$1 ":" \$2 "|" \$3; count[key]=\$4 } END { print "circRNA_ID", "count" > "${meta.id}.bsj.tsv"; for(i=1;i<=n;i++) { k=ids[i]; print k, (k in count ? count[k] : 0) >> "${meta.id}.bsj.tsv" } }' ${circ_list} ${result}.BSJ_Matrix
    awk 'BEGIN { FS=OFS="\\t" } NR==FNR { key=\$1 ":" \$2 "|" \$3; ids[++n]=key; next } FNR>1 && NF>=4 { key=\$1 ":" \$2 "|" \$3; count[key]=\$4 } END { print "circRNA_ID", "count" > "${meta.id}.fsj.tsv"; for(i=1;i<=n;i++) { k=ids[i]; print k, (k in count ? count[k] : 0) >> "${meta.id}.fsj.tsv" } }' ${circ_list} ${result}.FSJ_Matrix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>&1 | sed -n '1s/.*Awk //; s/,.*//p')
    END_VERSIONS
    """

    stub:
    """
    printf 'circRNA_ID\\tcount\\n' > ${meta.id}.bsj.tsv
    printf 'circRNA_ID\\tcount\\n' > ${meta.id}.fsj.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
