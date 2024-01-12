#!/usr/bin/perl
use strict; use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
BEGIN { unshift @INC, dirname(abs_path($0))."/scripts"; }
#FAlite: http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm
use FAlite;

#use DataBrowser;

die"
Usage: create_custom_transcriptome_fa.pl <reference genome fasta file> <annotation gtf file> <custom_transcriptome_fa output fasta file name>
" unless @ARGV==3;

my ($ref_fa, $annot_gtf, $output_file) = ($ARGV[0], $ARGV[1], $ARGV[2]);

for my $fa_file ($ref_fa, $output_file) {
	die "Error with $fa_file, fasta file names need to end with .fa or .fasta\n" unless ($fa_file =~ /.*(\.fa|\.fasta)$/);
}

my %transcripts_data;
my %chrs;

incr_transcripts_data($annot_gtf);

my %chr_fa;

incr_genome_fa();

open my $OUT, ">$output_file" or die "can't open $output_file\n";

for my $transcript_name (keys %transcripts_data) {
	my @exons;
#	if ($transcripts_data{$transcript_name}{strand} eq "-") {
#		@exons = sort {$b <=> $a} keys %{$transcripts_data{$transcript_name}{exons}};
#	} else {
#		@exons = sort {$a <=> $b} keys %{$transcripts_data{$transcript_name}{exons}};
#	}
	my %tmp;
	while (my ($key, $value) = each (%{$transcripts_data{$transcript_name}{exons}})){
		my @tmp_info = split(" ", $value);
		$tmp{$key}=$tmp_info[0];			
	}
	@exons = sort {$tmp{$a} <=> $tmp{$b}} keys %tmp;
#	@exons = sort {$transcripts_data{$transcript_name}{exons}{$a}[0] <=> $transcripts_data{$transcript_name}{exons}{$b}[0]} keys %{$transcripts_data{$transcript_name}{exons}}; #added by Christina to replace the above code block, in order to avoid bugs when using non-genecode annotations.
	
	my $first_exon = $exons[0];
	my @first_exon_info = split(" ", $transcripts_data{$transcript_name}{exons}{$first_exon});
	my $transcript_start = $first_exon_info[0];
	
	my $last_exon = $exons[-1];
	my @last_exon_info = split(" ", $transcripts_data{$transcript_name}{exons}{$last_exon});
	my $transcript_end = $last_exon_info[1];
	
	my @exon_lengths;
	my @exon_offs;
	for my $exon (@exons) {
		my $n_exon = $transcripts_data{$transcript_name}{exons}{$exon};
		my @n_exon_info = split(" ", $n_exon);
		my ($exon_start,$exon_end) = ($n_exon_info[0],$n_exon_info[1]);
		my $length = $exon_end - $exon_start + 1;
		my $off = $exon_start - $transcript_start;
		push @exon_lengths, $length;
		push @exon_offs, $off;
	}
	
	$transcript_start = $transcript_start-1;
	
	my $Txsequence;
	for (my $i=0; $i < @exon_lengths; $i++) {
		my $ExStart_pos = $transcript_start + $exon_offs[$i];
		$Txsequence .= substr($chr_fa{$transcripts_data{$transcript_name}{chr}},$ExStart_pos,$exon_lengths[$i]);
	}
	
	print $OUT ">" . $transcripts_data{$transcript_name}{chr} . ":$transcript_start-$transcript_end" . ":" . $transcripts_data{$transcript_name}{strand} .
	"\n" . $Txsequence . "\n";
}

close $OUT;

sub incr_genome_fa {

	open IN, "<$ref_fa" or die "can't open $ref_fa\n";
	my $fasta = new FAlite(\*IN);
	 while(my $entry = $fasta->nextEntry) {
	  	my $header = $entry->def;
	   	my $seq = $entry->seq;
	   	my ($ref_sequence_name) = $header =~ /^>(.*)/;
   	
	   	$chr_fa{$ref_sequence_name} .= uc($seq);
	}
	close IN;
}

sub incr_transcripts_data {
	my $gtf_file = shift;
	open my $IN, "<$gtf_file" or die "can't open $gtf_file\n";
	while (<$IN>) {
		if ($_ =~ /^#/) {
			next;
		}
		my @line = split("\t", $_);
		my ($chr,$element,$start,$end,$strand,$info) = ($line[0],$line[2],$line[3],$line[4],$line[6],$line[8]);
		next if ($element ne "exon");
		
		my ($gene_name) = $info =~ /gene_name "(.*?)";/;
		my ($transcript_id) = $info =~ /transcript_id "(.*?)";/;
		my ($exon_num) = $info =~ /exon_number (\d+)/;

		# Set the exon number to 1 if it is not defined
		$exon_num = 1 if (!$exon_num);
		
		my $transcript_name;
		if (!$gene_name) {
			$transcript_name = $transcript_id;
		} else {
			$transcript_name = $transcript_id . ";" . $gene_name;
		}
		
		$transcripts_data{$transcript_name}{chr} = $chr;
		$transcripts_data{$transcript_name}{strand} = $strand;
		$transcripts_data{$transcript_name}{exons}{$exon_num} = $start . " " . $end;
		
		$chrs{$chr}++;		
	}
	close $IN;
	
}