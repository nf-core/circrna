include { CIRCTOOLS_LONGREAD } from '../longread_tools/circtools'

workflow LONGREAD {
    take:
    ch_reads
    genome

    main:

    CIRCTOOLS_LONGREAD(ch_reads, genome)

    emit:
    reference = CIRCTOOLS_LONGREAD.out.reference
}
