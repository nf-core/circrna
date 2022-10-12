#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/circrna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/circrna

    Website: https://nf-co.re/circrna
    Slack  : https://nfcore.slack.com/channels/circrna
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CIRCRNA } from './workflows/circrna'

//
// WORKFLOW: Run main nf-core/circrna analysis pipeline
//
workflow NFCORE_CIRCRNA {
    CIRCRNA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CIRCRNA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MAPSPLICE_PARSE{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}/*",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/MapSplice/intermediates/${it}" : null }

    when:
    'mapsplice' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(raw_fusion) from mapsplice_fusion
    file(fasta) from ch_fasta
    file(gene_annotation) from ch_gene

    output:
    tuple val(base), file("${base}_mapsplice.bed") into mapsplice_results
    tuple val(base), file("${base}/*") into mapsplice_intermediates
    tuple val(base), val("MapSplice"), file("${base}_mapsplice_circs.bed") into mapsplice_annotated

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t MapSplice $raw_fusion -b ${base}/${base}.mapsplice.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}/${base}.mapsplice.junction.bed -o ${base}/${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_mapsplice.bed

    ## Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_mapsplice.bed > ${base}_mapsplice_circs.bed
    """
}

process SEGEMEHL_ALIGN{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/Segemehl/intermediates/${it}" : null }

    when:
    'segemehl' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(fastq) from segemehl_reads
    file(fasta) from ch_fasta
    file(idx) from ch_segemehl

    output:
    tuple val(base), file("${base}") into segemehl_intermediates
    tuple val(base), file("${base}_segemehl.bed") into segemehl_results
    tuple val(base), val("Segemehl"), file("${base}_segemehl_circs.bed") into segemehl_annotated

    script:
    def handleSam = params.save_quantification_intermediates ? "samtools view -hbS ${base}/${base}.sam > ${base}/${base}.bam && rm ${base}/${base}.sam" : "rm -rf ${base}/${base}.sam"
    """
    mkdir -p ${base}

    segemehl.x \\
        -t ${task.cpus} \\
        -d $fasta \\
        -i $idx \\
        -q ${fastq[0]} \\
        -p ${fastq[1]} \\
        -S \\
        -o ${base}/${base}.sam

    $handleSam

    ## Segemehl does not preserve strand information, nor account for it
    ## when collapsing and counting reads using haarz.x. This is my own fix which does.
    grep ';C;' ${base}/${base}.sngl.bed | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6}' | sort | uniq -c | awk -v OFS="\t" '{print \$2,\$3,\$4,\$5,\$1}' > ${base}/${base}_collapsed.bed

    ## Let user filter by BSJ read count param.
    awk -v OFS="\t" -v BSJ=${params.bsj_reads} '{if(\$5>=BSJ) print \$0}' ${base}/${base}_collapsed.bed > ${base}_segemehl.bed

    ## Re-work for Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_segemehl.bed > ${base}_segemehl_circs.bed
    """
}


// collect circrna info from each alignment process
ch_circs = ciriquant_annotated.mix(circexplorer2_annotated, dcc_annotated, circrna_finder_annotated, find_circ_annotated, mapsplice_annotated, segemehl_annotated)


process ANNOTATION{
    tag "${base}:${tool}"
    label 'process_high'
    publishDir "${params.outdir}/circrna_discovery/${tool}/${base}", mode: params.publish_dir_mode, pattern: "${base}.bed"
    publishDir "${params.outdir}/circrna_discovery/${tool}/annotation_logs", mode: params.publish_dir_mode, pattern: "${base}.log"

    input:
    tuple val(base), val(tool), file(circs) from ch_circs
    file(gtf_filt) from ch_gtf_filtered

    output:
    tuple val(base), val(tool), file("${base}.bed") into ch_annotation, ch_mirna_targets_bed12
    tuple val(base), file("${base}.log") into annotation_logs

    script:
    """
    mv $circs circs.bed
    bash ${workflow.projectDir}/bin/annotate_outputs.sh &> ${base}.log
    mv master_bed12.bed ${base}.bed.tmp

    ## isolate exon blocks
    awk -v FS="\t" '{print \$11}' ${base}.bed.tmp > mature_len.tmp

    ## sum exon block values
    awk -v FS="," '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' mature_len.tmp > mature_length

    ## concat to annotation file.
    paste ${base}.bed.tmp mature_length > ${base}.bed
    """

}


process FASTA{
    tag "${base}:${tool}"
    publishDir "${params.outdir}/circrna_discovery/${tool}/${base}", mode: params.publish_dir_mode, pattern: "*.fasta"

    input:
    tuple val(base), val(tool), file(bed) from ch_annotation
    file(fasta) from ch_fasta

    output:
    tuple val(base), val(tool), file("${base}.fa") into ch_mature_len_miranda, ch_mature_len_targetscan
    file("${base}.fasta") into publish_circ_mature_seq

    script:
    """
    ## FASTA sequences (bedtools does not like the extra annotation info - split will not work properly)
    cut -d\$'\t' -f1-12 ${base}.bed > bed12.tmp
    bedtools getfasta -fi $fasta -bed bed12.tmp -s -split -name > circ_seq.tmp
    ## clean fasta header
    grep -A 1 '>' circ_seq.tmp | cut -d: -f1,2,3 > ${base}.fa && rm circ_seq.tmp
    ## add backsplice sequence for miRanda TargetScan, publish canonical FASTA to results.
    rm $fasta
    bash ${workflow.projectDir}/bin/backsplice_gen.sh ${base}.fa
    """
}


/*
================================================================================
                        Generate circRNA count matrix
================================================================================
*/

quantification_results = ciriquant_results.mix(circexplorer2_results, circrna_finder_results, dcc_results, find_circ_results, mapsplice_results, segemehl_results)

if(tools_selected > 1){
    process MERGE_TOOLS{
        tag "${base}"

        when:
        'circrna_discovery' in module

        input:
        tuple val(base), file(bed) from quantification_results.groupTuple()

        output:
        file("${base}.bed") into sample_counts

        script:
        """
        ## make list of files for R to read
        ls *.bed > samples.csv

        ## Add catch for empty bed file and delete
        bash ${workflow.projectDir}/bin/check_empty.sh

        ## Use intersection of "n" (params.tool_filter) circRNAs called by tools
        ## remove duplicate IDs, keep highest count.
        Rscript ${workflow.projectDir}/bin/consolidate_algorithms_intersection.R samples.csv $params.tool_filter

        mv combined_counts.bed ${base}.bed
        """
    }
    process COUNT_MATRIX_COMBINED{
        publishDir "${params.outdir}/circrna_discovery", pattern: "count_matrix.txt", mode: params.publish_dir_mode

        when:
        'circrna_discovery' in module

        input:
        file(bed) from sample_counts.collect()

        output:
        file("circRNA_matrix.txt") into circRNA_counts
        file("count_matrix.txt") into matrix

        script:
        """
        python ${workflow.projectDir}/bin/circRNA_counts_matrix.py > matrix.txt
        ## handle non-canon chromosomes here (https://stackoverflow.com/questions/71479919/joining-columns-based-on-number-of-fields)
        n_samps=\$(ls *.bed | wc -l)
        canon=\$(awk -v a="\$n_samps" 'BEGIN {print a + 4}')
        awk -v n="\$canon" '{ for (i = 2; i <= NF - n + 1; ++i) { \$1 = \$1"-"\$i; \$i=""; } } 1' matrix.txt | awk -v OFS="\t" '\$1=\$1' > circRNA_matrix.txt
        Rscript ${workflow.projectDir}/bin/reformat_count_matrix.R
        """
    }
}else{
    process COUNT_MATRIX_SINGLE{
        publishDir "${params.outdir}/circrna_discovery", pattern: "count_matrix.txt", mode: params.publish_dir_mode

        when:
        'circrna_discovery' in module

        input:
        file(bed) from quantification_results.collect()
        val(tool) from params.tool

        output:
        file("circRNA_matrix.txt") into circRNA_counts
        file("count_matrix.txt") into matrix

        script:
        """
        # Strip tool name from BED files (no consolidation prior to this step for 1 tool)
        for b in *.bed; do
            basename=\${b%".bed"};
            sample_name=\${basename%"_${tool}"};
            mv \$b \${sample_name}.bed
        done

        python ${workflow.projectDir}/bin/circRNA_counts_matrix.py > matrix.txt
        ## handle non-canon chromosomes here (https://stackoverflow.com/questions/71479919/joining-columns-based-on-number-of-fields)
        n_samps=\$(ls *.bed | wc -l)
        canon=\$(awk -v a="\$n_samps" 'BEGIN {print a + 4}')
        awk -v n="\$canon" '{ for (i = 2; i <= NF - n + 1; ++i) { \$1 = \$1"-"\$i; \$i=""; } } 1' matrix.txt | awk -v OFS="\t" '\$1=\$1' > circRNA_matrix.txt
        Rscript ${workflow.projectDir}/bin/reformat_count_matrix.R
        """
    }
}

/*
================================================================================
                            miRNA Prediction
================================================================================
*/

process TARGETSCAN_DATABASE{
    when:
    'mirna_prediction' in module

    input:
    file(mature) from ch_mature

    output:
    file("mature.txt") into ch_mature_txt

    script:
    """
    bash ${workflow.projectDir}/bin/targetscan_format.sh $mature
    """
}


process MIRANDA{
    tag "${base}"
    label 'process_long'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.miRanda.txt",
        saveAs: { params.save_mirna_predictions ? "mirna_prediction/miRanda/${tool}/${base}/${it}" : null }
    when:
    'mirna_prediction' in module

    input:
    tuple val(base), val(tool), file(fasta) from ch_mature_len_miranda
    file(mirbase) from ch_mature
    file(mirbase_txt) from ch_mature_txt

    output:
    tuple val(base), val(tool), file("*.miRanda.txt") into ch_miranda_results

    script:
    prefix = fasta.toString() - ~/.fa/
    """
    miranda $mirbase $fasta -strict -out ${prefix}.bindsites.out -quiet
    echo "miRNA Target Score Energy_KcalMol Query_Start Query_End Subject_Start Subject_End Aln_len Subject_Identity Query_Identity" | tr ' ' '\t' > ${prefix}.miRanda.txt

    # Add catch here for non hits (supply NAs to outfile)
    # Making the decision that if miRanda fails, then the miRNA analysis for this circRNA exits cleanly.
    # Happy to rework in the future, but do not want pipeline failing on circrnas with no MRE sites.
    # targetscan produces far more results, have not observed it producing empty reults.
    ## exit code 1 = fail, 0 = success
    if grep -A 1 -q "Scores for this hit:" ${prefix}.bindsites.out;
    then
        grep -A 1 "Scores for this hit:" ${prefix}.bindsites.out | sort | grep ">" | cut -c 2- | tr ' ' '\t' >> ${prefix}.miRanda.txt
    else
        ## Add NA's to miRanda cols:
        printf "%0.sNA\t" {1..11} >> ${prefix}.miRanda.txt
    fi
    """
}


process TARGETSCAN{
    tag "${base}"
    label 'process_long'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.targetscan.txt",
        saveAs: { params.save_mirna_predictions ? "mirna_prediction/TargetScan/${tool}/${base}/${it}" : null }
    when:
    'mirna_prediction' in module

    input:
    tuple val(base), val(tool), file(fasta) from ch_mature_len_targetscan
    file(mirbase_txt) from ch_mature_txt

    output:
    tuple val(base), val(tool), file("*.targetscan.txt") into ch_targetscan_results

    script:
    prefix = fasta.toString() - ~/.fa/
    """
    ##format for targetscan
    cat $fasta | grep ">" | sed 's/>//g' > id
    cat $fasta | grep -v ">" > seq
    paste id seq | awk -v OFS="\t" '{print \$1, "0000", \$2}' > ${prefix}_ts.txt

    # run targetscan
    targetscan_70.pl mature.txt ${prefix}_ts.txt ${prefix}.targetscan.txt
    """
}

ch_mirna_targets = ch_miranda_results.join(ch_targetscan_results, by:[0,1])

ch_targets = ch_mirna_targets.join(ch_mirna_targets_bed12, by:[0,1])

process MIRNA_TARGETS{
    tag "${base}"
    label 'process_low'
    publishDir "${params.outdir}/mirna_prediction/${tool}", mode: params.publish_dir_mode, pattern: "*miRNA_targets.txt"

    input:
    tuple val(base), val(tool), file(miranda), file(targetscan), file(bed12) from ch_targets

    output:
    tuple val(base), file("*miRNA_targets.txt") into circrna_mirna_targets

    script:
    """
    ## reformat and sort miRanda, TargetScan outputs, convert to BED for overlaps.
    tail -n +2 $targetscan | sort -k1,1 -k4n | awk -v OFS="\t" '{print \$1, \$2, \$4, \$5, \$9}' | awk -v OFS="\t" '{print \$2, \$3, \$4, \$1, "0", \$5}' > targetscan.bed
    tail -n +2 $miranda | sort -k2,2 -k7n | awk -v OFS="\t" '{print \$2, \$1, \$3, \$4, \$7, \$8}' | awk -v OFS="\t" '{print \$2, \$5, \$6, \$1, \$3, \$4}' | sed 's/^[^-]*-//g' > miranda.bed

    ## intersect, consolidate miRanda, TargetScan information about miRs.
    ## -wa to output miRanda hits - targetscan makes it difficult to resolve duplicate miRNAs at MRE sites.
    bedtools intersect -a miranda.bed -b targetscan.bed -wa > ${base}.mirnas.tmp
    bedtools intersect -a targetscan.bed -b miranda.bed | awk '{print \$6}' > mirna_type

    ## remove duplicate miRNA entries at MRE sites.
    ## strategy: sory by circs, sort by start position, sort by site type - the goal is to take the best site type (i.e rank site type found at MRE site).
    paste ${base}.mirnas.tmp mirna_type | sort -k3,3 -k2n -k7r | awk -v OFS="\t" '{print \$4,\$1,\$2,\$3,\$5,\$6,\$7}' | awk -F "\t" '{if (!seen[\$1,\$2,\$3,\$4,\$5,\$6]++)print}' | sort -k1,1 -k3n > ${base}.mirna_targets.tmp
    echo -e "circRNA\tmiRNA\tStart\tEnd\tScore\tEnergy_KcalMol\tSite_type" | cat - ${base}.mirna_targets.tmp > ${base}.miRNA_targets.txt
    """
}

/*
================================================================================
                            Differential Expression
================================================================================
*/

// converseley to CIRIquant using Hisat2, this proc does need the index files. stage as new channel
ch_hisat_index_files = params.hisat ? Channel.value(file("${params.hisat}/*")) :  hisat_built

process HISAT_ALIGN{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}.bam",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/Hisat2/${it}" : null }

    when:
    'differential_expression' in module

    input:
    tuple val(base), file(fastq) from hisat_reads
    file(hisat2_index) from ch_hisat_index_files.collect()
    file(fasta) from ch_fasta

    output:
    tuple val(base), file("${base}.bam") into hisat_bam

    script:
    """
    hisat2 -p ${task.cpus} --dta -q -x ${fasta.baseName} -1 ${fastq[0]} -2 ${fastq[1]} -t | samtools view -bS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam
    """
}

process STRINGTIE{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/StringTie/${it}" : null }

    when:
    'differential_expression' in module

    input:
    tuple val(base), file(bam) from hisat_bam
    file(gtf) from ch_gtf

    output:
    file("${base}") into stringtie_dir

    script:
    """
    mkdir ${base}/
    stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p ${task.cpus} -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
    """
}

process DEA{
    label 'process_medium'
    publishDir "${params.outdir}/differential_expression", pattern: "circRNA", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "boxplots", mode: params.publish_dir_mode
    publishDir "${params.outdir}/quality_control", pattern: "DESeq2_QC", mode: params.publish_dir_mode
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "RNA-Seq",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/${it}" : null }

    when:
    'differential_expression' in module

    input:
    file(gtf_dir) from stringtie_dir.collect()
    file(circ_matrix) from circRNA_counts
    file(phenotype) from ch_phenotype
    val(species) from ch_species

    output:
    file("RNA-Seq") into rnaseq_dir
    file("circRNA") into circrna_dir
    file("boxplots") into boxplots_dir
    file("DESeq2_QC") into qc_plots

    script:
    """
    for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done

    prepDE.py -i samples.txt

    ## prepDE && circRNA counts headers are sorted where uppercase preceedes lowercase i.e Z before a
    ## reformat the phenotype file to match the order of the samples.
    head -n 1 $phenotype > header
    tail -n +2 $phenotype | LC_COLLATE=C sort > sorted_pheno
    cat header sorted_pheno > tmp && rm phenotype.csv && mv tmp phenotype.csv

    Rscript ${workflow.projectDir}/bin/DEA.R gene_count_matrix.csv $phenotype $circ_matrix $species ${workflow.projectDir}/bin/ensemblDatabase_map.txt

    mv gene_count_matrix.csv RNA-Seq
    mv transcript_count_matrix.csv RNA-Seq
    """
}

/*
================================================================================
                                    MultiQC
================================================================================
*/

process MULTIQC{
    label 'py3'
    label 'process_low'
    publishDir "${params.outdir}/quality_control/MultiQC", mode: params.publish_dir_mode,
        pattern: "*.html"

    input:
    file(raw_fastqc) from fastqc_raw.collect().ifEmpty([])
    file(trim_fastqc) from fastqc_trimmed.collect().ifEmpty([])
    file(BBDUK_stats) from bbduk_stats_ch.collect().ifEmpty([])
    file(multiqc_config) from ch_multiqc_config
    file(mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file("*.html") into multiqc_out

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}

/*
================================================================================
                            Auxiliary functions
================================================================================
*/

// Check integer
def isValidInteger(value){
    value instanceof Integer
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Define list of available tools
def defineToolList() {
    return [
        'ciriquant',
        'circexplorer2',
        'find_circ',
        'circrna_finder',
        'dcc',
        'mapsplice',
        'segemehl'
        ]
}

// Define module list
def defineModuleList() {
    return [
    'circrna_discovery',
    'mirna_prediction',
    'differential_expression'
    ]
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core/circrna] error:  Invalid CSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core/circrna] error: Cannot find supplied FASTQ or BAM input file. If using input method CSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Read input files from input CSV
def extract_data(csvFile){
    Channel
        .fromPath(csvFile)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_keys = ["Sample_ID", "Read1", "Read2", "Bam"]
        if(!row.keySet().containsAll(expected_keys)) exit 1, "[nf-core/circrna] error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2', 'Bam'."

        checkNumberOfItem(row, 4)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)
        def bam = row.Bam.matches('NA') ? 'NA' : return_file(row.Bam)

        if(samples == '' || read1 == '' || read2 == '' || bam == '') exit 1, "[nf-core/circrna] error: a field does not contain any information. Please check your CSV file"
        if(read1.matches('NA') && read2.matches('NA') && bam.matches('NA')) exit 1, "[nf-core/circrna] error: A row in your CSV file appears to have missing information."
        if( !read1.matches('NA') && !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "[nf-core/circrna] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
        if( !read2.matches('NA') && !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "[nf-core/circrna] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
        if( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified BAM file either has a non-recognizable extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

        // output tuple mimicking fromFilePairs if FASTQ provided, else tuple for BAM
        if(bam.matches('NA')){
            [ samples, [read1, read2] ]
        }else{
            [ samples, bam ]
        }

        }
}

// If no input CSV provided, parse input directory containing files.
def retrieve_input_paths(input, type){

    if(type == 'fastq'){

        fastq_files = input
        Channel
                .fromFilePairs(fastq_files)
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty{exit 1, "[nf-core/circrna] error: Your FASTQ files do not have the appropriate extension of either '.fastq.gz', '.fq.gz', .fastq' or '.fq'."}
                .map{ row -> [ row[0], [ row[1][0], row[1][1] ]]}
                .ifEmpty{exit 1, "[nf-core/circrna] error: --input was empty - no files supplied"}
                .set{reads_for_csv}

    }else if(type == 'bam'){

        bam_files = input
        Channel
                .fromFilePairs(bam_files, size: 1)
                .filter{ it =~/.*.bam/}
                .map{ row -> [row[0], [row[1][0]]]}
                .ifEmpty{exit 1, "[nf-core/circrna] error: Cannot find bam file matching: ${bam_files}"}
                .set{reads_for_csv}
    }

    reads_for_csv
                .map{

                def samples = it[0]
                def read1 = (type == 'bam') ? 'NA' : return_file(it[1][0])
                def read2 = (type == 'bam') ? 'NA' : return_file(it[1][1])
                def bam =   (type == 'fastq') ? 'NA' : return_file(it[1][0])

                if(bam.matches('NA')){
                    [ samples, [read1, read2] ]
                }else{
                    [ samples, bam ]
                }

                }
                .ifEmpty{exit 1, "[nf-core/circrna] error: Invalid file paths with --input"}
}


// Check input phenotype file

def examine_phenotype(pheno){

    Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['condition']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "[nf-core/circrna] error: 'condition' is not a column name in the phenotype file.\n\nThe response variable must be named 'condition', please refer to the usage documentation online"

        def condition  = row.condition.matches('NA') ? 'NA' : row.condition

        if(condition == '') exit 1, "[nf-core/circrna] error: Invalid phenotype file, condition column contains empty cells."
        if(condition.matches('NA')) exit 1, "[nf-core/circrna] error: NA value in phenotype condition column."

        }
        .toList()

        return Channel.value(file(pheno))
}

/*
================================================================================
                            nf-core functions
================================================================================
*/

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/circrna v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

// handle multiqc_channels
//if(params.trim_fastq == false){
//    ch_multiqc_report = multiqc_trim_out
//}else{
//    ch_multiqc_report = multiqc_raw_out
//}

// Completion e-mail notification
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/circrna] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/circrna] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/circrna] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/circrna] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$workflow.projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$workflow.projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$workflow.projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$workflow.projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
            mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/circrna]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/circrna]${c_red} Pipeline completed with errors${c_reset}-"
    }
}



def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
