#!/usr/bin/perl

use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/sequence_helpers.pl";

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

my $key_name = get_arg("k", "EXTRACT_ALL_KEYS", \%args);
my $start_location = get_arg("s", 1, \%args);
my $end_location = get_arg("e", -1, \%args);
my $cur_strand = get_arg("st", "", \%args);
my $locations_file = get_arg("f", "", \%args);
my $print_location_description = get_arg("d", 0, \%args);
my $print_location_name = get_arg("dn", 0, \%args);
my $print_location_description_chr = get_arg("dc", 0, \%args);
my $zero_based_coordinates = get_arg("0", 0, \%args);
my $gap = get_arg("g", 1, \%args);
my $minplus = get_arg("mp", 0, \%args);

if ($cur_strand ne "") { $minplus = 1; }

if ($gap<1) { $gap=1 }

my $subtraction_coordinates = $zero_based_coordinates == 1 ? 0 : 1;

my %key_starts;
my %key_ends;
my %key_strands;
my %key_output_names;
if (length($locations_file) > 0)
{
    open(LOCATIONS, "<$locations_file");
    while(<LOCATIONS>)
    {
	chop;

	my @row = split(/\t/);

	$key_starts{$row[0]} .= "$row[2]\t";
	$key_ends{$row[0]} .= "$row[3]\t";
	$key_output_names{$row[0]} .= "$row[1]\t";
	if ($minplus)
	{
		$key_strands{$row[0]} .= "$row[4]\t";
	}
	
	#print STDERR "key_starts{$row[0]}=$key_starts{$row[0]}\n";
	#print STDERR "key_ends{$row[0]}=$key_ends{$row[0]}\n";
	#print STDERR "key_output_names{$row[0]}=$key_output_names{$row[0]}\n";
    }
}
else
{
    $key_starts{$key_name} = $start_location;
    $key_ends{$key_name} = $end_location;
    $key_strands{$key_name} = $cur_strand;
    $key_output_names{$key_name} = $key_name;
}

while(<$file_ref>)
{
    chop;

    my @row = split(/\t/);

    if ((length($locations_file) == 0 and $key_name eq "EXTRACT_ALL_KEYS") or length($key_starts{$row[0]}) > 0)
    {
	my $key = length($key_starts{$row[0]}) > 0 ? $row[0] : "EXTRACT_ALL_KEYS";
	my @starts = split(/\t/, $key_starts{$key});
	my @ends = split(/\t/, $key_ends{$key});
	my @strands = split(/\t/, $key_strands{$key});
	my @output_names = split(/\t/, $key_output_names{$key});

	for (my $i = 0; $i < @starts; $i++)
	{
	    my $start_location = $starts[$i];
	    my $end_location = $ends[$i];
	    my $strand = $strands[$i];

	    my $output_name = length($output_names[$i]) > 0 ? $output_names[$i] : $row[0];

	    if ($print_location_description == 1)
	    {

		print "$output_name from $start_location to $end_location\t";
	    }
	    elsif ($print_location_name == 1)
	    {
		print "$output_name\t";
	    }
	    elsif ($print_location_description_chr == 1){
	      print "$row[0]\t$output_name\t$start_location\t$end_location\t";
	    }
	    else
	    {
		print "$row[0]\t";
	    }

	    if ($start_location == 0 and $zero_based_coordinates != 1) { $end_location = length($row[1]); }
	    if ($end_location == -1) { $end_location = length($row[1]); }

	    my $sequence_start = $start_location < $end_location ? $start_location : $end_location;

	    my $string = substr($row[1], $sequence_start - $subtraction_coordinates, abs($end_location - $start_location) + 1);

	    if ($gap>1) {
		my $new_string;
		for (my $j=0;$j<length($string);$j++){
		    if ($j-(int($j/$gap)*$gap) < 1) { $new_string .= substr($string,$j,1) }
		}
		$string=$new_string;
	    }

	    if ($minplus)
	    {
	    	if ($strand eq "-")
	    	{
	    		$string = &ReverseComplement($string);
	    	}
	    }
	    else
	    {
	    	if ($start_location > $end_location)
			{
				$string = &ReverseComplement($string);
			}
		}
	    print "$string\n";
	}
    }
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/extract_sequence.pl <file>

   Extracts a sequence from a given stab file

   NOTE: Locations here are specified in 1-based coordinates (change with -0)

   -0:       Zero-based coordinates (default: 1-based coordinates)

   -k <num>:  Key to extract (default: extract from all keys)
   -s <num>:  Start location (default: 1)
   -e <num>:  End location (default: end of sequence)
   -st <str>: Strand ('"+"' or '"-"') to force, otherwise determined by start>end.

   -f <str>:  File containing the locations to extract in the format:
              key<tab>name<tab>start<tab>end

   -mp:       Take the fifth column of the input location files as the strand (+ or -)
              and disregard the ordering of start, end. Good especially when extracting
              a sequence of length 1 (strand is anbigous when start=end).

   -d:        Print the full description of the location of the extracted sequence
   -dn:       Print only the name of the extracted sequence
   -dc:       Print description in chr format

   -g <num>:  Extract only every nth letter (default: 1).

