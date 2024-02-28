process PSIRC_COMBINE {
    label 'process_single'

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    path(files)

    output:
    path("counts.tsv"), emit: est_counts
    path("tpm.tsv"), emit: tpm


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

    counts_df.to_csv("counts.tsv", sep="\\t")
    tpm_df.to_csv("tpm.tsv", sep="\\t")
    """
}
