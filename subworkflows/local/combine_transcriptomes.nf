include { GNU_SORT as COMBINE_TRANSCRIPTOME_GTFS } from '../../modules/nf-core/gnu/sort'
include { GAWK as EXCLUDE_OVERLONG_TRANSCRIPTS   } from '../../modules/nf-core/gawk'
include { GFFREAD as TRANSCRIPTOME               } from '../../modules/nf-core/gffread'

workflow COMBINE_TRANSCRIPTOMES {
    take:
    ch_genome_fasta
    ch_genome_gtf
    ch_circ_gtf

    main:
    ch_versions = channel.empty()

    COMBINE_TRANSCRIPTOME_GTFS(
        ch_genome_gtf.mix(ch_circ_gtf).map{ _meta, gtf -> gtf }.collect().map{ gtfs -> [[id: "transcriptome"], gtfs]},
    )

    EXCLUDE_OVERLONG_TRANSCRIPTS(
        COMBINE_TRANSCRIPTOME_GTFS.out.sorted, [], false
    )

    TRANSCRIPTOME(
        EXCLUDE_OVERLONG_TRANSCRIPTS.out.output,
        ch_genome_fasta.map{_meta, fasta -> fasta}
    )

    emit:
    fasta = TRANSCRIPTOME.out.gffread_fasta
    gtf   = EXCLUDE_OVERLONG_TRANSCRIPTS.out.output

    versions = ch_versions
}
