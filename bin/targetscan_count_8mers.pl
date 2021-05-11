#! /usr/bin/env perl

# Find locations of TargetScan seeds in a sequence file
# Also print file of sequence lengths

if (! $ARGV[1])
{
	print STDERR "\nFind defined sites from a set of UTRs\n";
	print STDERR "USAGE: $0 miRNA_seeds_file UTRs > Offset-6mer_sites\n";	
	print STDERR "Note: Include last argument of '-list' to print list of all positions too\n\n";
	exit;
}

###  Example miR-1 : seed region = GGAAUGU
###  To find and count all 6mer-1a (offset 6mer) sites,
###  miRNA family should be listed as uGGAAU 
###  (initial T + positions 2-6 of the mature miRNA)
###  Then script will look for the reverse complement

$miRNASeedsFile = $ARGV[0];
$sequenceFile = $ARGV[1];

if ($sequenceFile =~ /(\S+)\.txt$/)
{
	$sequenceLengthsOutfile = "$1.lengths.txt";
}
else
{
	$sequenceLengthsOutfile = "$sequenceFile.lengths.txt";
}

# Create output file to hold sequence lengths
open (SEQUENCE_LENGTHS_OUT, ">$sequenceLengthsOutfile") || die "Cannot open $sequenceLengthsOutfile for writing: $!";

open (DEFINED_SITES, "$miRNASeedsFile") || die "Cannot open $miRNASeedsFile for reading: $!";
while (<DEFINED_SITES>)
{
	# GGAAUGU	GGAAUGU	10090	# miR-1
	# AAAGCCA	AAAGCCA	10090
	# AAAGCCA	AAAGCCA	10116

	chomp;
	($familyID, $familyID, $species) = split /\t/, $_;	
	$seedRegionRevcomp = reverse_complement($familyID);
	# Convert any Ts to Us
	$seedRegionRevcomp =~ s/T/U/g;
	$seedRegionRevcomp =~ s/t/u/g;

	$siteToFind = $seedRegionRevcomp . "A";
	
	$familyID2site{$familyID} = $siteToFind;
	
	# Reverse-complement to get the actual site we need
	push @sitesToFind, $siteToFind;
}

# Make these unique, since multiple families can have the same 6mer-1a site
@sitesToFind = do { my %seen; grep { !$seen{$_}++ } @sitesToFind };

open (SEQUENCES, "$sequenceFile") || die "Cannot open $sequenceFile for reading: $!";
while (<SEQUENCES>)
{
	# ENST00000000233	9713	-CAACCAGGGGCCG------G-CCCCTGCTGC
	
	chomp;
	($sequenceID, $species, $sequence) = split /\t/, $_;
	
	# Remove all gaps
	$sequence =~ s/-//g;
	$sequence =~ s/\.//g;
	# Convert any Ts to Us
	$sequence =~ s/T/U/g;
	$sequence =~ s/t/u/g;
	
	$sequenceLength = length $sequence;
	print SEQUENCE_LENGTHS_OUT "$sequenceID\t$species\t$sequenceLength\n";
	
	%siteToPositions = ();
	
	foreach $siteToFind (@sitesToFind)
	{
		$numThisSite = 0;
		
		# Make sure this is case-insensitive (i)
		while ($sequence =~ /$siteToFind/gi)
		{
			$start = $-[0];
			$end = $+[0];
			
			$numThisSite++;
		}
		$siteToCount{$sequenceID}{$species}{$siteToFind} = $numThisSite;
	}
	
	foreach $familyID (sort keys %familyID2site)
	{
		$siteToFind = $familyID2site{$familyID};
		if ($siteToCount{$sequenceID}{$species}{$siteToFind})
		{
			print "$sequenceID\t$species\t$familyID\t$siteToCount{$sequenceID}{$species}{$siteToFind}\n";
		}
	}
}

print STDERR "\nAll done -- Also see $sequenceLengthsOutfile for sequence lengths.\n\n";

##########################################

sub reverse_complement 
{
	my $dna = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;
	return $revcomp;
}
