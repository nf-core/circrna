process PSIRC_COMBINE {
    label 'process_single'

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    path(files)

    output:
    path("linear_counts.tsv"),   emit: linear_counts
    path("circular_counts.tsv"), emit: circular_counts
    path("counts.tsv"),          emit: counts

    path("linear_tpm.tsv"),      emit: linear_tpm
    path("circular_tpm.tsv"),    emit: circular_tpm
    path("tpm.tsv"),             emit: tpm


    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd

    files = "${files}".split(" ")

    counts_df = pd.DataFrame()
    tpm_df = pd.DataFrame()

    for f in files:
        sample_name = os.path.basename(f)[:-len(".tsv")]
        df = pd.read_csv(f, sep="\\t", index_col=0)

        counts_df[sample_name] = df["est_counts"]
        tpm_df[sample_name] = df["tpm"]

    for name, df in [("counts", counts_df), ("tpm", tpm_df)]:
        is_linear = df.index.str.startswith("ENS")
        linear_df = df[is_linear]
        circular_df = df[~is_linear]

        linear_df.to_csv(f"linear_{name}.tsv", sep="\\t")
        circular_df.to_csv(f"circular_{name}.tsv", sep="\\t")
        df.to_csv(f"{name}.tsv", sep="\\t")
    """
}
