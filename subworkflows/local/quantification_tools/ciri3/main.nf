include { CIRI3_CIRI3 as CIRI3 } from '../../../../modules/local/ciri3/ciri3'
include { CIRI3_UC_LIST } from '../../../../modules/local/ciri3/uc_list'
include { CIRI3_UC_MATRIX } from '../../../../modules/local/ciri3/uc_matrix'
include { JOIN_SAMPLES as JOIN_BSJ } from '../../../../modules/local/matrix/join_samples'
include { JOIN_SAMPLES as JOIN_FSJ } from '../../../../modules/local/matrix/join_samples'

workflow CIRI3_QUANT {
    take:
    sam
    bed12
    ch_fasta
    ch_gtf

    main:

    CIRI3_UC_LIST(bed12)

    ch_ciri_list = CIRI3_UC_LIST.out.list.first()
    ch_fasta_value = ch_fasta.first()
    ch_gtf_value = ch_gtf.first()
    CIRI3(sam, ch_fasta_value, ch_gtf_value, ch_ciri_list.map { _meta, list -> list })

    // Preserve sample association; CIRI3's txt result is optional in -UC mode.
    ch_matrix_pair = CIRI3.out.bsj_matrix.join(CIRI3.out.fsj_matrix)
    CIRI3_UC_MATRIX(ch_ciri_list, ch_matrix_pair)

    JOIN_BSJ(
        CIRI3_UC_MATRIX.out.bsj.map { meta, table -> [[id: 'ciri3_bsj'], meta.id, table] }.groupTuple(),
        'circRNA_ID'
    )

    JOIN_FSJ(
        CIRI3_UC_MATRIX.out.fsj.map { meta, table -> [[id: 'ciri3_fsj'], meta.id, table] }.groupTuple(),
        'circRNA_ID'
    )

    emit:
    bsj      = JOIN_BSJ.out.joined.map { _meta, table -> [[id: 'ciri3_bsj', quantification: 'ciri3'], table] }
    fsj      = JOIN_FSJ.out.joined.map { _meta, table -> [[id: 'ciri3_fsj', quantification: 'ciri3'], table] }
}
