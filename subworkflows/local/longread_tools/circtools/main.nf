include { CIRCTOOLS_REFERENCE } from '../../../../modules/local/circtools/reference'
include { CIRCTOOLS_NANOPORE } from '../../../../modules/local/circtools/nanopore'

workflow CIRCTOOLS_LONGREAD {
    take:
    ch_reads
    genome

    main:

    CIRCTOOLS_REFERENCE(genome)

    CIRCTOOLS_NANOPORE(ch_reads, CIRCTOOLS_REFERENCE.out.reference, genome)

    emit:
    reference = CIRCTOOLS_REFERENCE.out.reference
}
