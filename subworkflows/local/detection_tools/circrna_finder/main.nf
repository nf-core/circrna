include { CIRCRNA_FINDER as MAIN } from '../../../../modules/local/circrna_finder'
include { GAWK as UNIFY          } from '../../../../modules/nf-core/gawk'

workflow CIRCRNA_FINDER {
    take:
    star_sam
    star_junctions
    star_tab

    main:

    ch_joined = star_sam
        .join(star_junctions)
        .join(star_tab)
        .map { meta, sam, junction, tab ->
            [meta + [tool: "circrna_finder"], [sam, junction, tab]]
        }

    MAIN(ch_joined)

    UNIFY(MAIN.out.results.map { meta, results -> [meta, results, "circrna_finder.bed"] }, [], false)

    emit:
    bed = UNIFY.out.output
}
