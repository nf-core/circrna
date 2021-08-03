#!/usr/bin/env nextflow

/*
 * Below has been removed from main.nf (after diff_exp process)
 * reworked DEA.R to produce correct dirs + boxplots
 * Saved for posterity and tests after commit
 */

 (circrna_dir_fetch, circrna_dir_plots) = circrna_dir.into(2)

 // obtain the IDs of the differentially expressed circrna

 process fetch_de_circ_id{

        input:
            file(circrna_dir) from circrna_dir_fetch

        output:
            file("*.bed") into de_bed_files

        when: 'differential_expression' in module

        script:
        up_reg = "${circrna_dir}/*up_regulated_differential_expression.txt"
            down_reg = "${circrna_dir}/*down_regulated_differential_expression.txt"
        """
        grep -v "baseMean" $up_reg > up_reg_noheader.txt
        cat $down_reg up_reg_noheader.txt > de_circs.txt

        # make dummy files out of these to place in channel
        awk '{print \$1}' de_circs.txt | grep -v "ID" | while read -r line; do touch \${line}.bed12.bed; done
        """
 }

 // de circrna dummy bed files in de_bed_files channel, map to get simpleName
 ch_de_bed = de_bed_files.flatten().map{ file -> [file.simpleName, file]}

 // match incoming bed file channel of all circrnas with 'dummy' DE circ bed files
 filt_bed_tmp = ch_de_bed.join(bed_diff_exp)

 // remove the dummy files (file[1])
 filt_bed = filt_bed_tmp.map{file -> [file[0], file[2]]}

 // join should keep only the DE keys and apply it to parent, mature files
 ch_report = filt_bed.join(parent_diff_exp).join(mature_diff_exp)

 // must combine folders here or else process uses once then exits.
 ch_DESeq2_dirs = circrna_dir_plots.combine(rnaseq_dir)

 process de_plots{

        publishDir "${params.outdir}/differential_expression/circrna_expression_plots", pattern:"*.pdf", mode:'copy'

            input:
                file(phenotype) from ch_phenotype
                tuple val(base), file(bed), file(parent_gene), file(mature_length), file(circRNA), file(rnaseq) from ch_report.combine(ch_DESeq2_dirs)

            output:
                file("*.pdf") into de_plots
            file("*DESeq2_stats.txt") into de_stats

        when: 'differential_expression' in module

            script:
            up_reg = "${circRNA}/*up_regulated_differential_expression.txt"
            down_reg = "${circRNA}/*down_regulated_differential_expression.txt"
            circ_counts = "${circRNA}/DESeq2_normalized_counts.txt"
            gene_counts = "${rnaseq}/DESeq2_normalized_counts.txt"
            """
            # merge upreg, downreg info
        grep -v "baseMean" $up_reg > up_reg_noheader.txt
        cat $down_reg up_reg_noheader.txt > de_circ.txt

            # Make plots and generate circRNA info
            Rscript ${projectDir}/bin/circ_report.R de_circ.txt $circ_counts $gene_counts $parent_gene $bed $mature_length $phenotype
            """
 }

 // collect all from previous process
 //master_ch = circRNA_plots.collect()
 //(test, test1) = circRNA_plots.into(2)
 //test.view()
 // delete text files in process script, left with only dirs.

 process master_report{

        publishDir "${params.outdir}/differential_expression/circrna_diff_exp_stats", mode:'copy'

            input:
                file(reports) from de_stats.collect()

            output:
                file("*circRNAs.txt") into final_out

        when: 'differential_expression' in module

            script:
            """
            # remove header, add manually
            cat *.txt > merged.txt
            grep -v "Log2FC" merged.txt > no_headers.txt
            echo "circRNA_ID Type Mature_Length Parent_Gene Strand Log2FC pvalue Adjusted_pvalue" | tr ' ' '\t' > headers.txt
            cat headers.txt no_headers.txt > merged_reports.txt

            Rscript ${projectDir}/bin/annotate_report.R
            """
 }
