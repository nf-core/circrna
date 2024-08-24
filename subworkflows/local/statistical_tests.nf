include { CIRCTEST_PREPARE  } from '../../modules/local/circtest/prepare'
include { CIRCTEST_CIRCTEST } from '../../modules/local/circtest/circtest'
include { CIRIQUANT_PREPDE     } from '../../modules/local/ciriquant/prepde'
include { STRINGTIE_PREPDE  } from '../../modules/local/stringtie/prepde'

workflow STATISTICAL_TESTS {
    take:
    ch_gene_counts
    ch_circ_counts
    ch_ciriquant
    ch_stringtie
    ch_phenotype

    main:
    ch_versions = Channel.empty()

    ch_counts = ch_gene_counts.join(ch_circ_counts)

    CIRCTEST_PREPARE(ch_counts)
    ch_versions = ch_versions.mix(CIRCTEST_PREPARE.out.versions)

    CIRCTEST_CIRCTEST(CIRCTEST_PREPARE.out.counts,
        ch_phenotype)
    ch_versions = ch_versions.mix(CIRCTEST_CIRCTEST.out.versions)

    ch_phenotype_annotations = ch_phenotype
        .map{ meta, table -> table.text }
        .splitCsv( header: true )

    ch_condition_samples = ch_phenotype_annotations
        .map{ annotations -> [annotations.sample, annotations.condition] }
        .join(ch_ciriquant.map{ meta, gtf -> [meta.id, gtf] })
        .join(ch_stringtie.map{ meta, gtf -> [meta.id, gtf] })
        .groupTuple(by: 1)
        .map{ samples, condition, ciriquant, stringtie -> [condition, samples, ciriquant, stringtie] }

    ch_condition_pairs = ch_condition_samples
        .combine(ch_condition_samples)
        .filter{ c_control, s_control, f_ciri_control, f_stringtie_control, c_treatment, s_treatment, f_ciri_treatment, f_stringtie_treatment
             -> c_control != c_treatment }
        .map{ c_control, s_control, f_ciri_control, f_stringtie_control, c_treatment, s_treatment, f_ciri_treatment, f_stringtie_treatment ->
            [   [id: "${c_control}:${c_treatment}"],
                s_control + s_treatment,
                f_ciri_control + f_ciri_treatment,
                f_stringtie_control + f_stringtie_treatment,
                ['C'] * s_control.size() + ['T'] * s_treatment.size()]}

    CIRIQUANT_PREPDE(ch_condition_pairs
        .map{meta, samples, ciri, stringtie, conditions -> [meta, samples, ciri, conditions]}
    )
    // ch_versions = ch_versions.mix(CIRIQUANT_DEA.out.versions)

    STRINGTIE_PREPDE(ch_condition_pairs
        .map{meta, samples, ciri, stringtie, conditions -> [meta, samples, stringtie]}
    )
    //ch_versions = ch_versions.mix(STRINGTIE_PREPDE.out.versions)

    emit:
    versions = ch_versions
}
