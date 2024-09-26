#!/usr/bin/perl

use strict;
use List::Util qw[shuffle min max];

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/sequence_helpers.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/format_number.pl";

my $RNAHYBRID_EXE_DIR = "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/Bin/RNAHybrid/RNAhybrid-2.1/src";
my $RNAddG_EXE_DIR 		= "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/Bin/ViennaRNA/ViennaRNA-1.6/Progs/";


if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}

my $file_ref;
my $file = $ARGV[0];
if (length($file) < 1 or $file =~ /^-/) 
{
  $file_ref = \*STDIN;
}
else
{
  open(FILE, $file) or die("Could not open file '$file'.\n");
  $file_ref = \*FILE;
}

my %args = load_args(\@ARGV);

my $ddG_area = get_arg("ddgarea", "70", \%args);
my $FULL_TL  = get_arg("tl", "50", \%args);
my $DDG_OPEN  = get_arg("dgtl", "25", \%args);
my $include_3max = get_arg("3max", 0, \%args);
my $include_ddG_v4 = get_arg("ddG_v4", 0, \%args);
my $components = get_arg("components", "", \%args);

my $upstream_rest = get_arg("upstream_rest", 0, \%args);
my $downstream_rest = get_arg("downstream_rest", 0, \%args);
my $no_force = get_arg("noforce", 0, \%args);

# Step over all Locations ###########################################################

my $done = 0;
my $bunch_size = 1000;

while (!$done)
{
	my @dGduplexesInputArray = ();
	my @ddGInputArray        = ();
	my @headerInputArray     = ();
	my $record;
	
	for (my $nlines=0; $nlines < $bunch_size; $nlines++)
	{
		if (!($record = <$file_ref>))
		{
			$done = 1;
			last;
		}
	
		chomp ($record);
		my @row = split(/\t/, $record);
		my $endpos = $row[2];
		my $miRNA = $row[8];
		my $utr = $row[11];
		my $force_binding_i ="";
		my $force_binding_j ="";
	
		if (!$no_force)
		{
			# $force_binding = "-i $row[9] -j $row[10]";
			$force_binding_i = $row[9];
			$force_binding_j = $row[10];
		}
		
		# Calculate begining of target but make sure we're still in the UTR
		
		my $startpos  = max ($endpos - $FULL_TL + 1, 1);
		
		my $maxtargetLen = $endpos - $startpos + 1;
		my $target        = substr ($utr, $startpos - 1, $maxtargetLen);
		my $real_start = max ($endpos - $DDG_OPEN + 1, 1);
	
		push (@headerInputArray,     join ("\t", splice(@row, 0, 8)));
		push (@dGduplexesInputArray, join ("\t", $miRNA, $target, $force_binding_i, $force_binding_j));
		push (@ddGInputArray,        join ("\t", $miRNA, $utr, $real_start, $endpos, $upstream_rest, $downstream_rest, $ddG_area));
	}
	
	my  $arraySize = @headerInputArray;

	print STDERR "\nComputing $arraySize results: ";
	
	my @dGduplexesOutputArray = &getdGduplexes (@dGduplexesInputArray);
	my @ddGOutputArray = &getddG (@ddGInputArray);
		
	for (my $i=0; $i < $arraySize; $i++)
	{
			my ($target_length, $dGall, $dG5, $dG3) = split (/\t/, $dGduplexesOutputArray[$i]);
			my ($ddGall, $dG0, $dG1) = split (/\t/, $ddGOutputArray[$i]);

			$ddGall = $dGall + ($dG1 - $dG0);
			
			my $ddG_v4   = 0;		
			my $dG3max   = 0;
			my $dG3ratio = 0;
	
			print $headerInputArray[$i] . "\t$target_length\t$dGall\t$dG5\t$dG3\t$ddGall\t$dG0\t$dG1\t$ddG_v4\t$dG3max\t$dG3ratio\n";
	}
}

print STDERR " Done.\n";
# system ("rm tmp_seqfile1 tmp_seqfile2");


################################################################################

sub getdGduplexes {

	my @inArray = @_;
	my @outArray = ();
	
	open (SEQFILE, ">tmp_seqfile1") or die ("Could not open temporary sequence file.\n");
	
	my $insize = @inArray;
	
	foreach my $oneTarget (@inArray)
	{
		(my $miRNA, my $target, my $force_binding_i, my $force_binding_j) = split (/\t/, $oneTarget);
		print SEQFILE "$miRNA\n$target\n$force_binding_i\n$force_binding_j\n";
	}
	
	close (SEQFILE);


	# Call RNAduplex and extract result dGs and length
	
	my $cmd = "$RNAddG_EXE_DIR/RNAduplex -5 0 < tmp_seqfile1";
	print STDERR "Calling RNAduplex with " . @inArray . " targets... ";

	my $result_of_cmd = `$cmd`;
	
	my @resArray = split (/\n/, $result_of_cmd);
	
	my $outsize = @resArray;
	
	if ($insize ne $outsize)
	{
		die "RNAduplex failure. Result was $result_of_cmd\n";
	}
	
	foreach my $resline (@resArray)
	{
		my ($ret_structure, $ret_start_miR, $ret_end_miR, $ret_start_target, $ret_end_target, $ret_miR_len, $ret_target_len, $ret_dGall, $ret_dG5, $ret_dG3) = split (/\t/, $resline);
		push (@outArray, join ("\t", $ret_target_len, $ret_dGall, $ret_dG5, $ret_dG3));
	}
	
	return @outArray;
}


################################################################################

sub getddG {

	my @inArray = @_;
	my @outArray = ();
	my $upstream_rest;
	my @seq_area;

	open (SEQFILE, ">tmp_seqfile2") or die ("Could not open temporary sequence file.\n");
	
	foreach my $oneTarget (@inArray)
	{
		(my $miRNA, my $utr, my $startpos, my $endpos, $upstream_rest, my $downstream_rest, my $ddG_area) = split (/\t/, $oneTarget);

		# Extract area around sequence
		@seq_area = &extract_sequence ($startpos, $endpos, $downstream_rest + $ddG_area, $upstream_rest + $ddG_area, $utr);

		print SEQFILE "$seq_area[2]\n";
	}
	
	close (SEQFILE);

	if ($upstream_rest != "0")
	{
		$upstream_rest =~ s/^0*//;		# Because we're gonna call C with this argument and sscanf is problematic with leading zeros
	}
	
	my $seq = $seq_area[2];
	my $bindStart = $seq_area[0];
	my $bindEnd   = $seq_area[1];

	my $target_len = $bindEnd - $bindStart + 1 + $downstream_rest;
	
	my $cmd = "$RNAddG_EXE_DIR/RNAddG4 -u $upstream_rest -s $bindEnd -f $target_len -t $target_len < tmp_seqfile2";

	print STDERR "Calling RNAddG with " . @inArray . " targets... ";

	my @resArray = split (/\n/, `$cmd`);

	foreach my $resline (@resArray)
	{
		my @ddG_values = split (/\t/, $resline);
		push (@outArray, join ("\t", $ddG_values[6], $ddG_values[3], $ddG_values[4]));
	}
	
	return @outArray;
}

################################################################################

sub extract_sequence {

	my $start_location = $_[0];
	my $end_location = $_[1];
	my $cur_us = $_[2];
	my $cur_ds = $_[3];
	my $sequence = $_[4];

	#QQQ print STDERR "start = $start_location, end=$end_location; cur_us=$cur_us; cur_ds=$cur_ds\n";
	
	my $us_start = $start_location > $cur_us ? $start_location - $cur_us : 1;
	my $ds_end  = $end_location + $cur_ds < length($sequence) ? $end_location + $cur_ds : length($sequence);
	my $actual_start = $start_location - $us_start + 1;

	my $string = substr($sequence, $us_start - 1, ($ds_end - $us_start) + 1);
	
	my $actual_end = $actual_start + abs($end_location - $start_location);

	my @result = ($actual_start, $actual_end, $string);
	
	#QQQ print STDERR "string=$string length=" . length($string) . "\n";
	
	return @result;
}

# 
# ################################################################################
# 
# sub getdG {
# 
# 	print "getdG OBSOLETE. Please use getdGduplex.\n"; exit;
# 	
# 	my ($miRNA, $miR_start, $miR_len, $target, $target_start, $target_len) = @_ ;
# 
# 	# Build miR and target substrings
# 
# 	my $submiR    = substr ($miRNA, $miR_start-1, $miR_len);
# 	my $subtarget = reverse (substr (reverse($target), $target_start-1, $target_len)); 
# 
# 	# Call RNAHybrid and extract result dG
# 	
# 	my $prog_line = "$RNAHYBRID_EXE_DIR/RNAhybrid -c -s 3utr_fly -b 1 $subtarget $submiR";
# 	my @prog_result = split (/\n/, `$prog_line`);
# 	if (@prog_result != 1)
# 	{
# 		print STDERR "Cound not hybrid.\n";
# 		print STDERR "Program call was: $prog_line\n";
# 		print STDERR "Result was: " . join ("\n", @prog_result);
# 		exit;
# 	}
# 	
# 	
# #QQQ	print STDERR "Program call was: $prog_line\n";
# 	
# 	my @result = split (/:/, $prog_result[0]);
# 	my $dG = $result[4];
# 	
# 	return $dG;
# }
# 
# ################################################################################
# 
# sub getdGduplex {
# 
# 	print "getdGduplex OBSOLETE. Please use getdGduplexes.\n"; exit;
# 
# 	my ($miRNA, $miR_start, $miR_len, $target, $target_start, $target_len) = @_ ;
# 
# 	# Build miR and target substrings
# 
# 	my $submiR    = substr ($miRNA, $miR_start-1, $miR_len);
# 	my $subtarget = reverse (substr (reverse($target), $target_start-1, $target_len)); 
# 
# 	# Call RNAduplex and extract result dG
# 
# 	open (SEQFILE, ">tmp_seqfile1") or die ("Could not open temporary sequence file.\n");
# 	print SEQFILE "$subtarget\n$submiR\n";
# 	close (SEQFILE);
# 
# 	print STDERR "\nCalling duplex with subtarget=$subtarget; submiR=$submiR\n";
# 	
# 	my $resline = `$RNAddG_EXE_DIR/RNAduplex < tmp_seqfile1`;
# 	chomp ($resline);
# 	
# 	print STDERR "Result is: $resline\n";
# 
# 	$resline =~ m/ \((.*)\)/;
# 	return $1;
# }


__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/RNAddG_compute.pl <file>

	Compute the dG energies of potential target sites whose seed is given in the external
	tab file.
	
	Output consists of the input ID, followed by the following:
	    
        Type         miRNA coordinates           Target coordinates
        ====         =================           ==================
        target_length
        dG all       nt 1 and up                 nt 1-FULL_TL
        dG 5'        nt 1-9                      computed
        dG 3'        nt 10 and up                computed
        dG0
        dG1
        ddG_v4
        dG 3' max    nt 10 and up                Reverse complement of
                                                 miRNA bases 10 and up
        dG 3' ratio
		
    -tl <num>:           Maximal target length (default: 50)
    -dgtl <num>:         Target length when opening for ddG calculation (default: 25)
    -ddgarea <num>:      Area upstream and downstream around target to fold (default: 70) 
    -no3max:             Do not compute 3' max value and ratio (saves time)
    -components <file>:  Print the components of the ddG calculation into the given file.
    
    
 
