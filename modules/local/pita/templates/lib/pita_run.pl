#!/usr/bin/perl

use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";

if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}


my $CDS_FLANK = 200;
my $SITE_LENGTH = 25;

my $POLY_A = "A" x $CDS_FLANK;

my %args = load_args(\@ARGV);

my $prefix = get_arg("prefix", "", \%args);
my $debug = get_arg("debug", 0, \%args);
my $utr_fn = get_arg("utr", "", \%args);
my $upstream_fn = get_arg("upstream", "", \%args);
my $mir_fn = get_arg("mir", "", \%args);
my $FLANK_UP = get_arg("flank_up", 0, \%args);
my $FLANK_DOWN = get_arg("flank_down", 0, \%args);
my $output_dir = get_arg("output", "", \%args);
my $DDG_AREA = get_arg("ddG_context", 70, \%args);
my $use_stab = get_arg("stab", 0, \%args);
my $limit = get_arg("limit", "", \%args);

if ($output_dir ne "")
{
	$output_dir .= "/";
}


my $echoed_pt_args = 
  echo_arg_quoted("l", \%args) .
  echo_arg_quoted("gu", \%args) .
  echo_arg_quoted("m", \%args) .
  echo_arg_quoted("loop", \%args);
 
my $r = int(rand(100000));

if ($utr_fn eq "") { die "Must provide UTR fasta file.\n"; }
if ($mir_fn eq "") { die "Must provide microRNA fasta file.\n";}

## Step 0: Prepare extended UTR sequences

if ($use_stab)
{
	dsystem ("cp $utr_fn tmp_utr_stab_$r");
	dsystem ("cp $mir_fn tmp_mir_stab_$r");	
}
else
{
	dsystem ("cat $utr_fn | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/dos2unix.pl | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/fasta2stab.pl > tmp_utr_stab_$r");
	dsystem ("cat $mir_fn | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/dos2unix.pl | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/fasta2stab.pl > tmp_mir_stab_$r");
}

if ($upstream_fn eq "")
{
	dsystem ("cat tmp_utr_stab_$r " .
			 "| cut -f 1 " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/add_column.pl -s $POLY_A " .
			 "> tmp_prefix_$r;");
}
else
{
	if ($use_stab)
	{
		dsystem ("cp $upstream_fn tmp_upstream_stab_$r");
	}
	else
	{
		dsystem ("cat $upstream_fn | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/fasta2stab.pl > tmp_upstream_stab_$r");
	}
	
	dsystem ("cat tmp_utr_stab_$r " .
			 "| cut -f 1 " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/add_column.pl -s $POLY_A " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/join.pl -o A - tmp_upstream_stab_$r " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/merge_columns.pl -1 1 -2 2 " .
			 "> tmp_ext_utr_$r;");
	
	dsystem ("cat tmp_ext_utr_$r " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/stab2length.pl " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cut.pl -f 1,1,2,2 " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 2 -s " . ($CDS_FLANK-1) . " " .
			 ">tmp_ext_utr_pos_$r");
			 
	dsystem ("cat tmp_ext_utr_$r " .
			 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/extract_sequence.pl -f tmp_ext_utr_pos_$r " .
			 "> tmp_prefix_$r");
}

my $ext_utr_fn = $output_dir . $prefix . "ext_utr.stab";

dsystem ("/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/join.pl tmp_prefix_$r tmp_utr_stab_$r " .
         "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/add_column.pl -s $POLY_A " .
         "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/merge_columns.pl -1 1 -2 2 " .
         "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/merge_columns.pl -1 1 -2 2 " .
         "> " . $ext_utr_fn);



## Step 1: Search potential targets

print STDERR "Finding potential targets...\n";
dsystem ("/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/find_potential_mirna_targets.pl tmp_utr_stab_$r -f tmp_mir_stab_$r $echoed_pt_args > tmp_pt_$r");

my $n_potential_targets = `wc -l tmp_pt_$r | sed 's/^ *//g' | cut -f 1 -d ' '`;
if (length($limit) > 0 and $n_potential_targets > $limit)
{
   	dsystem ("/bin/rm -rf tmp_prefix_$r tmp_ext_utr_$r tmp_ext_utr_pos_$r tmp_prefix_$r tmp_pt_$r tmp_utr_stab_$r tmp_mir_stab_$r tmp_upstream_stab_$r");
	exit(1);
}

## Step 2: Compute site scores

print STDERR "Computing site scores...\n";
dsystem ("cat tmp_pt_$r " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/join.pl -1 2 - $ext_utr_fn " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cut.pl -f 1-9,11- " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 2,3 -a $CDS_FLANK " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/RNAddG_compute.pl -ddgarea $DDG_AREA -upstream_rest $FLANK_UP -downstream_rest $FLANK_DOWN " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 2,3 -s $CDS_FLANK " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 9,10,11,12,13,14 -m '\"-1\"' " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cut.pl -f 1-15,14 " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 15 -sc 14 " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cut.pl -f 1-12,14-16,13 " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl -c 9,10,11,12,13,14,15 -p 2 -m '\"-1\"' " .
		 "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cap.pl \"UTR,microRNA,Start,End,Seed,Mismatchs,G:U,Loop,Site size,dGduplex,dG5,dG3,dG0,dG1,dGopen,ddG\" " .
		 "> " . $output_dir . $prefix . "pita_results.tab");


## Clean up
if (not ($debug))
{
	dsystem ("/bin/rm -rf tmp_prefix_$r tmp_ext_utr_$r tmp_ext_utr_pos_$r tmp_prefix_$r tmp_pt_$r tmp_utr_stab_$r tmp_mir_stab_$r tmp_upstream_stab_$r");
}

exit (0);


sub dsystem {

        my $cmd = $_[0];

        if ($debug) { print STDERR "Executing $cmd\n"; }
        system ($cmd);
        
}


__DATA__
syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/pita_run.pl [OPTIONS]

Execute the PITA algorithm for identifying and scoring microRNA target sites.

options:
    -utr <filename>:      fasta file containing the UTRs to be scanned
    -mir <filename>:      fasta file containing the microRNA sequences
    -upstream <filename>: fasta file containing the upstream sequence for each UTR. The IDs
                          in should match the IDs found int the UTR file. If less 200 bp are
                          given (or if no file is given), it is padded with Poly-A.
    -stab:                All input files (utr, mir, upstream) are assumed to be stab files
                          rather than fasta files
                          
    -flank_up <bp>
    -flank_down <bp>:     Flank requirement in basepairs (default: zero for both)
    
    -ddG_context <bp>:    Number of bases upstream and downstream for target site that are
                          taken into account when folding the UTR (default: 70)
                          
    -output <directory>:  Path to where output files should be placed (default: current directory)
    -prefix <string>:     Add the string as a prefix to the output files (pita_results.tab and ext_utr.stab)
    -limit  <num>:        Limit the run to <num> potential targets (optional, no limit by default). Exits if limit is reached.
    -debug                Run in debug mode (leave untouched tmp files)

    
    
    Seed matching parameters:
    
    -l <num1-num2>:       Search for seed lengths of num1,...,num2 to the MicroRNA (default: 6-8)

    -gu <nums>:           Lengths for which G:U wobbles are allowed and number of allowed wobbles.
                          Format of nums: <length;num G:U>,<length;num G:U>,... (default: 6;0,7;1,8;1)

    -m <nums>:            Lengths for which mismatches are allowed and number of allowed mismatches
                          Format of nums: <length;num mismatches>,<length;num mismatches>,...
                          (default: 6;0,7;0,8;1)

    -loop <nums>:         Lengths for which a single loop in either the target or the microrna is allowed
                          Format of nums: <length>,<length>,... (default: none)
