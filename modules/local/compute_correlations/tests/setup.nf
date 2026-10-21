include { PSIRCQUANT_INDEX  } from '../../psircquant/index/main'
include { PSIRCQUANT_QUANT  } from '../../psircquant/quant/main'
include { TXIMETA_TXIMETA   } from '../../tximeta/tximeta/main'
include { MERGE_EXPERIMENTS } from '../../quantification/merge_experiments/main'

workflow CORRELATION_INPUT {

    main:
    def base = params.pipelines_testdata_base_path
    def transcriptome = channel.value([[id: 'transcriptome'], file(base + 'mini/intermediates/transcriptome.fasta', checkIfExists: true)])
    def gtf = channel.value([[id: 'transcriptome'], file(base + 'mini/intermediates/transcriptome.filtered.gtf', checkIfExists: true)])
    def fai = channel.value([[id: 'genome'], file(base + 'mini/intermediates/chrI_mini.fa.fai', checkIfExists: true)])
    def reads = channel.of(
            ['cel_N2_1', 'N2_rep1'],
            ['cel_N2_2', 'N2_rep2'],
            ['cel_N2_3', 'N2_rep3'],
            ['fust1_1', 'fust1_rep1'],
            ['fust1_2', 'fust1_rep2'],
            ['fust1_3', 'fust1_rep3'],
        )
        .map { sample, stem ->
            [[id: sample, single_end: false], [1, 2].collect { mate -> file(base + "fastq/${stem}_${mate}.fastq.gz", checkIfExists: true) }]
        }
    PSIRCQUANT_INDEX(transcriptome)
    PSIRCQUANT_QUANT(reads, PSIRCQUANT_INDEX.out.index, gtf, fai, 3)
    TXIMETA_TXIMETA(PSIRCQUANT_QUANT.out.directory, 'kallisto')
    // Small text table derived directly from quantification, without a test-only process.
    def tpm = PSIRCQUANT_QUANT.out.directory
        .collect(flat: false)
        .map { directories ->
            def samples = directories
                .sort { entry -> entry[0].id }
                .collectEntries { meta, directory ->
                    [(meta.id): directory.resolve('abundance.tsv').readLines().drop(1).collectEntries { line ->
                        def cells = line.split('\t')
                        [(cells[0]): cells[4]]
                    }]
                }
            def rows = [['transcript_id', 'gene_id'] + samples.keySet().toList()]
            samples
                .values()
                .first()
                .keySet()
                .sort()
                .each { transcript ->
                    rows << [transcript, transcript] + samples.values().collect { values -> values[transcript] }
                }
            rows.collect { row -> row.join('\t') }.join('\n') + '\n'
        }
        .collectFile(name: 'tpm.tsv')
        .map { table -> [[id: 'tpm'], table] }
    MERGE_EXPERIMENTS(
        TXIMETA_TXIMETA.out.se.map { _meta, se -> se }.collect().map { experiments -> [[id: 'six_samples'], experiments] },
        channel.value([[], []]),
        gtf,
        tpm,
    )

    emit:
    merged = MERGE_EXPERIMENTS.out.merged
}
