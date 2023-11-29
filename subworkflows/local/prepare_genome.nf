
include { BOWTIE_BUILD        } from '../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD       } from '../../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX           } from '../../modules/nf-core/bwa/index/main'
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
    ch_gtf   = Channel.fromPath(gtf)
    fasta_tuple = Channel.value([[id: "fasta"], fasta])
    gtf_tuple = Channel.value([[id: "gtf"], gtf])

    // MapSplice & find_circ requires reference genome to be split per chromosome:
    if( ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) && params.module.contains('circrna_discovery') ){
        directory = file("${params.outdir}/references/chromosomes")

        if ( !directory.exists() ) {
            directory.mkdirs()
            ch_fasta.splitFasta( record: [id:true] )
                .map{ record -> record.id.toString() }
                .set{ ID }

            ch_fasta.splitFasta( file: true )
                .merge( ID ).map{ file, id -> file.copyTo(directory + "/${id}.fa") }
        }

        stage_chromosomes = Channel.value(directory)
    }

    BOWTIE_BUILD(ch_fasta)

    BOWTIE2_BUILD(fasta_tuple)

    BWA_INDEX (fasta_tuple)

    HISAT2_EXTRACTSPLICESITES(gtf_tuple)

    HISAT2_BUILD(fasta_tuple, gtf_tuple, HISAT2_EXTRACTSPLICESITES.out.txt)

    STAR_GENOMEGENERATE(fasta_tuple, gtf_tuple)

    SEGEMEHL_INDEX(fasta)

    // Collect versions
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
    ch_versions = ch_versions.mix(SEGEMEHL_INDEX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    bowtie       = BOWTIE_BUILD.out.index
    bowtie2      = BOWTIE2_BUILD.out.index
    bwa          = BWA_INDEX.out.index
    chromosomes  = ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) ? stage_chromosomes : 'null'
    hisat2       = HISAT2_BUILD.out.index.collect()
    star         = STAR_GENOMEGENERATE.out.index.collect()
    segemehl     = SEGEMEHL_INDEX.out.index
    splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt.collect()
    versions     = ch_versions
}
