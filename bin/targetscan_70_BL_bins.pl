#! /usr/bin/env perl
use warnings;

#######################################################################
# Copyright(c) 2009-2015 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 7.0 (August 2015)
#
# Comment: This program uses set(s) of UTR alignments to get branch length (BL) at each position
#          and then the median across all positions.
#          This BL is used to assign UTR to a bin (1 - 10) based on degree of conservation.
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70
#
# Version 6 includes a new set of BL thresholds; otherwise, it's the same as version 6.
# Version 7 includes more species, a new tree, and a new set of BL thresholds; otherwise, it's the same as version 6.
#
#######################################################################

# UTR alignment input file (also used for TargetScan prediction code)
$UTRfile = $ARGV[0];	# Name of parsed MAF file

# Need this to read species tree
#   so BioPerl needs to be installed
use Bio::TreeIO;

# Need this to get median
#   so http://search.cpan.org/dist/Statistics-Lite/ needs to be installed
use Statistics::Lite qw(:all);

# Generic tree to use to get branch lengths for binning purposes
# $treeFile = "pct_trees_parameters/Unbinned_3UTR_tree.txt";
$treeFile = "PCT_parameters/Tree.generic.txt";

$refGenome = 9606;

# Minimum BL for each each bin (1 - 10) based on TargetScan 5.1 3' UTR partitioning into ten groups
# @BL_thresholds = qw (0 0.65655 0.97574 1.18258 1.35518 1.52027 1.68867 1.87412 2.11711 2.49238 );	# Release 5
# @BL_thresholds = qw (0 0.78088 1.15743 1.40192 1.59294 1.76514 1.93737 2.12238 2.37091 2.77733 );	# Release 5
@BL_thresholds = qw (0 1.21207417 2.17396073 2.80215414 3.26272822 3.65499277 4.01461968 4.40729032 4.90457274 5.78196252);

# What are all the nodes in the tree?
# These are NCBI Taxonomy IDs included in the analysis
@all_orgs = qw (7897 8364 8469 8478 8496 8839 8932 8954 9031 9103 9258 9305 9315 9361 9365 9371 9402 9483 9541 9544 9557 9595 9598 9601 9606 9615 9646 9669 9685 9708 9713 9733 9739 9785 9796 9807 9823 9913 9925 9940 9978 9986 10029 10036 10090 10116 10141 10160 10181 13146 13616 13735 27679 28377 28737 29078 30538 30611 34839 42254 43179 44394 48883 51337 55534 59463 59538 59729 59894 60711 61853 79684 127582 132908 143302 176014 181119 185453 225400 241585 246437 345164 419612 1230840);

$VERBOSE = "";

# Note: We can safely ignore the several errors from Bio/Tree/Node.pm at the end of the analysis.

###########################################  End of main constants  ###########################################

if (! $ARGV[0])
{
	print STDERR "\nCalculate median branch length from each UTR alignment\n";
	print STDERR "  and then assign to a UTR bin.\n";
	print STDERR "\nUSAGE: $0 tabbedAlignmentFile > Gene_BL_bin_file\n";
	print STDERR "Ex: $0 UTR_Sequences_sample.txt > UTRs_median_BLs_bins.txt\n";
	print STDERR "\ntabbedAlignmentFile fields: <seqID> <species> <sequence>\n";
	print STDERR "    where <species> belongs to the set of NCBI Taxonomy IDs\n";
	print STDERR "    (such as 9606 for human) used by TargetScanHuman and Mouse 5.\n";
	print STDERR "\nGene_BL_bin_file fields: <seqID> <branch_length> <bin[1-10]>\n\n";

	exit;
}

# Read the unbinned 3' UTR tree 
READ_TREE($treeFile);

# What thresholds are we using to assign a UTR to a bin (using BL)?
GET_BL_THRESHOLDS();

# Initialize
$lastGeneID = " ";

open (ALIGNMENT_TABS, $UTRfile) || die "Cannot open $UTRfile for reading: $!";
while (<ALIGNMENT_TABS>)
{
	# Read the UTR alignment

	# CDC2L6	9606	CCCACUCCCU---CU---------

	chomp;
	s/\r//g;
	
	@f = split (/\t/, $_);
	
	if ($geneID && $f[0] ne $lastGeneID)	# Process the last block of sequences
	{
		if ($VERBOSE)
		{
			print STDERR "Processing $lastGeneID...\n";
		}
		
		# Process an alignment for a UTR (gene) to get median BL and then bin
		getSpeciesListBLwithRefNt(0);
		
		$medianBLthisUTR = processUTRbranchLengths($lastGeneID);
		$binThisUTR = assignBLtoBin($medianBLthisUTR);
		
		print "$lastGeneID\t$medianBLthisUTR\t$binThisUTR\n";
		
		# Reset
		$geneID = "";
		@species = ();
		%speciesToFinalAlignment = ();
	}
	
	$geneID = $f[0];
	$species2alignment{$f[1]} = "$f[2]";
		
	push @species, $f[1];
	
	$lastGeneID = $geneID;
}

###  Get the last one
getSpeciesListBLwithRefNt(0);

$medianBLthisUTR = processUTRbranchLengths($lastGeneID);
$binThisUTR = assignBLtoBin($medianBLthisUTR);

print "$lastGeneID\t$medianBLthisUTR\t$binThisUTR\n";


#######################  Subroutines  #######################

sub getSpeciesListBLwithRefNt
{
	# Get a list of species with the same nt at this position as the reference genome

	my $printDetails = shift;

	$alignmentLength = length ($species2alignment{$refGenome});
	
	for ($i = 0; $i < $alignmentLength; $i++)
	{
		$keep = 0;
		$ntsThisPos = "";
	
		foreach $species (@species)
		{
			$posThisSpecies{$species} = substr($species2alignment{$species}, $i, 1);
			
			if ($posThisSpecies{$species} ne "-")
			{
				$keep++;
			}
		}
		
		if ($keep)
		{
			foreach $species (@species)
			{
				$speciesToFinalAlignment{$species} .= $posThisSpecies{$species};
				
				# Make a string of all nts at this position
				$ntsThisPos .= $posThisSpecies{$species};				
			}
						
			#
			#
			#
			# Set reference nt according to ref genome 
			#
			#
			#
			
			$consensusThisPos = "$posThisSpecies{$refGenome}"; 
			
			# Skip positions with a gap in the reference genome
			if  ($consensusThisPos ne "-")
			{
				# Compare each nt to consensus and count the ones that agree

				$consensusForBranchLength = uc($consensusThisPos);

				@speciesForBranchCalc = ();

				foreach $species (@species)
				{
					$ntThisPos = $posThisSpecies{$species};

					if ( lc($ntThisPos) eq lc($consensusForBranchLength))
					{
						push @speciesForBranchCalc, $species;
					}
				}
				
				# Sort all species in this list (to reduce number of combinations)
				@speciesForBranchCalc = sort @speciesForBranchCalc;
				
				$speciesThisPos = "@speciesForBranchCalc";

				###  Get branch length for this set of species at this position

				if ($speciesThisPos !~ /\s+/)	# Only one species in list
				{
					$branchLength = 0;
				}
				elsif ($speciesListToBL{$speciesThisPos})
				{
					$branchLength = $speciesListToBL{$speciesThisPos};
				}
				else
				{
					@include_orgs = split (/ /, $speciesThisPos);

					$branchLength = get_branch_length(@include_orgs);

					# Save in memory
					$speciesListToBL{$speciesThisPos} = $branchLength;
				}
				push @branchLengthsThisUTR, $branchLength;

				if ($printDetails)
				{
					print "$geneID\t$i\t$consensusForBranchLength\t$speciesThisPos\t$branchLength\n";
				}
			}
		}
	}
}

sub READ_TREE 
{
	# Read a precomputed phylogenetic tree

	my $treefile = shift;
	my $input = new Bio::TreeIO(-file   => $treefile,
					-format => "newick");
	$tree = $input->next_tree;

	foreach my $org (@all_orgs) 
	{
		($nodes{$org}) = $tree->find_node(-id => $org);
	}
}

sub GET_BL_THRESHOLDS
{	
	# Make s hash of BL thresholds
	
	$bin = 1;
	
	foreach $BL_threshold (@BL_thresholds)
	{
		$maxBLThisBin{$BL_threshold} = $bin;
		$bin++;
	}
}

sub get_branch_length 
{
	# Calculate branch length from a list of species IDs
	#
	# This section of code largely written by Robin Friedman, MIT

	my @orgs = @_;
	
	my %ref_ancestors;	# keys are internal ids, values are objects
	my %ref_cumul_dist;	# keys are internal ids, values 

	# are cumulative distance from node1 to given node
	my $place = $nodes{$refGenome};		# start at node1
	my $cumul_dist = 0;
	
	my @save = ();
	while ( $place )
	{
		my $id = $place->internal_id;
		push(@save, $id);

		@{$ref_ancestors{$id}} = @save;
		$ref_cumul_dist{$id} = $cumul_dist;
		if ($place->branch_length) 
		{
			$cumul_dist += $place->branch_length; # include current branch
		}
		
		$place = $place->ancestor;
	}
	
	# now climb up node2, for each node checking whether 
	# it's in node1_ancestors
	my %included_nodes = ();
	my $total_dist = 0;

	foreach my $org (@orgs) 
	{
		$place = $nodes{$org};	# start at leaf node
		$cumul_dist = 0;

		@save = ();
		while ( $place )
		{
			$id = $place->internal_id;
			push(@save, $id);

			if($included_nodes{$id}) 
			{
				$total_dist += $cumul_dist;
		
				last;
			}
			if(defined $ref_ancestors{$id}) 
			{
				# we're at lca
				$total_dist += $ref_cumul_dist{$id} + $cumul_dist;
	
				for (my $i = @{$ref_ancestors{$id}}-1; $i >= 0; $i--) 
				{
					my $cur = ${$ref_ancestors{$id}}[$i];

					if($included_nodes{$cur}) 
					{
						$total_dist -= $ref_cumul_dist{$cur};

						last;
					}
					else 
					{
						$included_nodes{$ref_ancestors{$id}[$i]} = 1;
					}
				}
				
				$included_nodes{$id} = 1;
				last;
			}
			$included_nodes{$id} = 1;

			# include current branch length in next iteration
			$cumul_dist += $place->branch_length || 0;
			
			$place = $place->ancestor;
		}
	}
	
	return $total_dist;
}

sub processUTRbranchLengths
{
	# Get the median of all branch lengths for a UTR

	my $medianBLthisUTR = median(@branchLengthsThisUTR);
	
	@branchLengthsThisUTR = ();
	
	return $medianBLthisUTR;
}

sub assignBLtoBin
{
	# Use BL thresholds to assign a UTR to a bin

	my $BL = shift;
	my $binthisBL;
	
	$binthisBL = 1;
	
	my @thresholds = sort { $a <=> $b } keys %maxBLThisBin;
	
	foreach my $threshold (@thresholds) 
	{
		if ($BL > $threshold)
		{
			$binthisBL = $maxBLThisBin{$threshold};
		}
	}
	
	return $binthisBL;
}
