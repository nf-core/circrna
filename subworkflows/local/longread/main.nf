include { CIRCTOOLS_LONGREAD } from '../longread_tools/circtools'

workflow LONGREAD {
    take:
    ch_reads
    genome

    main:

    CIRCTOOLS_LONGREAD(ch_reads, genome)

}
