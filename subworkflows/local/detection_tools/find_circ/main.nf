include { BOWTIE2_ALIGN as ALIGN       } from '../../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_VIEW                } from '../../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX               } from '../../../../modules/nf-core/samtools/index'
include { FIND_CIRC_ANCHORS as ANCHORS } from '../../../../modules/local/find_circ/anchors'
include { FIND_CIRC as MAIN            } from '../../../../modules/local/find_circ/find_circ'
include { BIOAWK as EXTRACT_READS      } from '../../../../modules/nf-core/bioawk'
include { CSVTK_SUMMARY as GROUP_READS } from '../../../../modules/local/csvtk/summary'
include { CSVTK_JOIN as JOIN_READS     } from '../../../../modules/nf-core/csvtk/join'
include { GAWK as UNIFY                } from '../../../../modules/nf-core/gawk'

workflow FIND_CIRC {
    take:
    reads
    bowtie2_index
    ch_fasta

    main:

    ALIGN( reads, bowtie2_index, ch_fasta, false, true )

    SAMTOOLS_INDEX( ALIGN.out.bam )

    SAMTOOLS_VIEW(
        ALIGN.out.bam.join( SAMTOOLS_INDEX.out.index ),
        ch_fasta.map { meta, fasta -> [meta, fasta, []] },
        [ [id: 'qname'], [] ],
        [ [id: 'bed'], [] ],
        ''
    )

    ANCHORS( SAMTOOLS_VIEW.out.bam )

    ch_find_circ = ANCHORS.out.anchors.multiMap { meta, anchors ->
        anchors:      [meta, anchors]
        strandedness: meta.strandedness ?: 'auto'
    }

    MAIN( ch_find_circ.anchors, bowtie2_index, ch_fasta.map{ _meta, fasta -> fasta }, ch_find_circ.strandedness )

    EXTRACT_READS( MAIN.out.reads, [], false, "fastq.gz" )

    GROUP_READS( EXTRACT_READS.out.output )

    JOIN_READS( MAIN.out.bed.join(GROUP_READS.out.csv)
        .map{ meta, bed, _reads -> [ meta, [bed, _reads]] }
    )

    UNIFY( JOIN_READS.out.out_file.map{ meta, bed ->
        [ meta + [tool: "find_circ"], bed, "find_circ.bed" ] }, [], false )

    emit:
    bed = UNIFY.out.output

}
