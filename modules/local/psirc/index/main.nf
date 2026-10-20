process PSIRC_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("*.all_possible_bsj_targets.fa.index"), path("*.FirstandLastExons_entities.fa.index"), emit: index
    tuple val("${task.process}"), val('psirc'), val('1.0.0'), emit: versions_psirc, topic: versions
    tuple val("${task.process}"), val('kallisto_b'), eval('kallisto_b version | sed "s/.*version //"'), emit: versions_kallisto_b, topic: versions

    script:
    """
    set +e
    psirc -i ${fasta} "\$(command -v kallisto_b)" > psirc-index.log 2>&1
    status=\$?
    set -e
    case \$status in
        0|255) ;;
        *) cat psirc-index.log >&2; exit \$status ;;
    esac
    if ! grep -q 'Finished generating' psirc-index.log || [ ! -s ${fasta.baseName}.all_possible_bsj_targets.fa.index ] || [ ! -s ${fasta.baseName}.FirstandLastExons_entities.fa.index ]; then
        cat psirc-index.log >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${fasta.baseName}.all_possible_bsj_targets.fa.index
    touch ${fasta.baseName}.FirstandLastExons_entities.fa.index
    """
}
