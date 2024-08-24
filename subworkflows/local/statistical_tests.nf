include { CIRCTEST_PREPARE  } from '../../modules/local/circtest/prepare'
include { CIRCTEST_CIRCTEST } from '../../modules/local/circtest/circtest'

workflow STATISTICAL_TESTS {
    take:
    ch_gene_counts
    ch_circ_counts
    ch_ciriquant
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
        .groupTuple(by: 1)
        .map{ samples, condition, files -> [condition, samples, files] }

    ch_condition_pairs = ch_condition_samples
        .combine(ch_condition_samples)
        .filter{ c_control, s_control, f_control, c_case, s_case, f_case -> c_control != c_case }
        .map{ c_control, s_control, f_control, c_case, s_case, f_case ->
            [s_control + s_case, f_control + f_case, [c_control] * s_control.size() + [c_case] * s_case.size()] }
        .view()

    emit:
    versions = ch_versions
}
