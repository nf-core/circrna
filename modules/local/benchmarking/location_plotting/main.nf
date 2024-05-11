process PLOT_LOCI {
    tag "Plotting loci from ${bed.baseName}"

    input:
    tuple val(meta), path(bed)

    output:
    path "*.png", emit: plots

    script:
    """
    python plot_bed.py ${bed} ${bed.baseName}.png
    """
}
