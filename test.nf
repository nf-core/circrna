#!/usr/bin/env nextflow

params.index = 'https://raw.githubusercontent.com/BarryDigby/test/dev/test_samples.csv'

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleID, [file(row.read1), file(row.read2)]) }
    .set{ samples_ch }


( foo_reads, raw_reads ) = samples_ch.into(2)


process foo{

	echo true

	input:
	tuple val(base), file(reads) from foo_reads

  shell:
  '''
  echo "+++++++++++++++++++"
  echo "Stats for !{base} read1"
  zcat !{reads[0]} | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
  echo "Stats for !{base} read2"
  zcat !{reads[1]} | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
  echo "+++++++++++++++++++"
  '''
}

aligner_reads = raw_reads

aligner_reads.view()
