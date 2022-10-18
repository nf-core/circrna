
include { BOWTIE_BUILD        } from '../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD       } from '../../modules/nf-core/bowtie2/build/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD        } from '../../modules/nf-core/hisat2/build/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { SEGEMEHL_INDEX      } from '../../modules/nf-core/segemehl/index/main'

workflow PREPARE_GENOME {

    take:
    fasta
    gtf

    main:
    ch_versions = Channel.empty()

    ch_fasta = Channel.fromPath(fasta)

    // MapSplice & find_circ requires reference genome to be split per chromosome:
     if( ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) && params.module.contains('circrna_discovery') ){
        file("${params.outdir}/genome/chromosomes").mkdirs()
        ch_fasta.splitFasta(record: [id:true])
                .map{ record -> record.id.toString() }
                .flatten()
                .set{ ID }

        ch_fasta.splitFasta(file: true)
                .flatten()
                .join(ID).map{ it ->
                             file = it[0]
                             chr_id = it[1]
                             file.copyTo("${params.outdir}/genome/chromosomes/${chr_id}.fa")
                       }

    stage_chromosomes = Channel.value("${workflow.launchDir}/${params.outdir}/genome/chromosomes")
    }

    // some index procs use tuple, some dont -_-
    ch_fasta.map{ it ->
             meta = [:]
             meta.id = it.simpleName
             return [ meta, [it] ]
    }.set{ fasta_tuple }

    BOWTIE_BUILD(
        fasta
    )

     BOWTIE2_BUILD(
        fasta_tuple
    )

    HISAT2_EXTRACTSPLICESITES(
        gtf
    )

    HISAT2_BUILD(
        fasta,
        gtf,
        HISAT2_EXTRACTSPLICESITES.out.txt
    )

    STAR_GENOMEGENERATE(
        fasta,
        gtf
    )

    SEGEMEHL_INDEX(
        fasta
    )

    // Collect versions
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(SEGEMEHL_INDEX.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first().ifEmpty(null))

    emit:
    bowtie      = BOWTIE_BUILD.out.index
    bowtie2     = BOWTIE2_BUILD.out.index
    chromosomes = ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) ? stage_chromosomes : 'null'
    hisat2      = HISAT2_BUILD.out.index
    star        = STAR_GENOMEGENERATE.out.index
    segemehl    = SEGEMEHL_INDEX.out.index
    versions    = ch_versions
}
