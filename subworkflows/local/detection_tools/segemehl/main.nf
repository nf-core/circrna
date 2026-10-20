include { SEGEMEHL_INDEX as INDEX   } from '../../../../modules/nf-core/segemehl/index'
include { SEGEMEHL_ALIGN as ALIGN   } from '../../../../modules/nf-core/segemehl/align'
include { GAWK as EXTRACT           } from '../../../../modules/nf-core/gawk'
include { GNU_SORT as SORT          } from '../../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_GROUPBY as GROUP } from '../../../../modules/nf-core/bedtools/groupby'
include { GAWK as UNIFY             } from '../../../../modules/nf-core/gawk'

workflow SEGEMEHL {
    take:
    reads
    fasta
    index

    main:

    // The workflow accepts an optional index path (as supplied by the pipeline),
    // while the upstream module consumes metadata-bearing tuples. Normalize both
    // reference inputs to tuples and pair by position: their metadata IDs need not match.
    if (fasta instanceof Path || fasta instanceof String) {
        ch_fasta = channel.value([[id: 'segemehl_reference'], file(fasta)])
    }
    else {
        ch_fasta = fasta.map { entry ->
            if (entry instanceof List) {
                entry
            }
            else {
                [[id: 'segemehl_reference'], entry]
            }
        }
    }

    if (index) {
        if (index instanceof Path || index instanceof String) {
            ch_index = channel.value([[id: 'segemehl_reference'], file(index)])
        }
        else {
            ch_index = index
        }
    }
    else {
        ch_index = INDEX(ch_fasta).index
    }

    ch_reference = ch_fasta
        .combine(ch_index)
        .map { fasta_meta, reference_fasta, _index_meta, reference_index -> [fasta_meta, reference_fasta, reference_index] }

    ALIGN(reads, ch_reference.first())
    EXTRACT(
        ALIGN.out.single_bed.map { meta, bed -> [meta + [tool: "segemehl"], bed, "segemehl_extracted.bed"] },
        [],
        false,
    )

    SORT(EXTRACT.out.output.map { meta, bed -> [meta, bed, "segemehl_sorted.bed"] })
    GROUP(SORT.out.sorted, 5)
    UNIFY(GROUP.out.bed.map { meta, bed -> [meta, bed, "segemehl.unified.bed"] }, [], false)

    emit:
    bed = UNIFY.out.output
}
