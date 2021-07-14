#! /usr/bin/env perl

#######################################################################
# Copyright(c) 2007-2015 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: George Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Version: 7.0
#
# Comment: This program predicts miRNA targets using the TargetScanS algorithm.
#          It produces output as displayed in TargetScan
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70
#
# Version: 1.1 corrects for 8mer 'UTR end' OBOB
# Version: 5.0 (11 Nov 2008)
#          - allows grouping of sites that are 8mer, 7mer-m8 AND 7mer-1a
#          - These are defined as a new type: "8mer+m8+1a"
#          - Group type is now a list of types in the group, rather than what's conserved
#          - An additional column in printed to show species in the group with a species' site type
#          - Handling of memory for large datasets has been improved
# Version: 5.1 -- Skipped, so code version will match TargetScan Release
# Version: 5.2 (21 July 2010)
#          - Corrected bug in method to identify 8mers that include 7mer-m8 and/or 7mer-1a sites
#          - When a site is found, modified backtrack so it's always far enough to get overlapping sites
# Version: 5.3/6.0 (Aug/Sept 2010)
#          - Add three more site types and overhaul some subroutines
# Version: 7.0 (July 2014 - June 2015)
#          - Skip target prediction in ribosome shadow region
#          - Look for sites only in miRNA-family-dependent species.
#          - Look for sites in ORF-overlapping regions (lowercase UTR) but flag as such
#
#######################################################################

# 28 Sep 2015 -- Modify reading of miRNA families to work with lowercase seed regions

# Basic ideas:
#
# 1 - Grab all miRNA info.
# 2 - Read through UTRs, getting those from one gene at a time.
# 3 - Identify miRNA sites for this gene.
# 4 - Group overlapping miRNA sites in different species into a group
#

use warnings;
use strict;	# Added in version 5.3

# Find sites in species even in which the miRNA has not been annotated
# If you don't want to identify these, set the variable to 0.  For TS7, this is set to 0.
our $FIND_SITES_ALL_SPECIES = 0;
# For site comparison between species, how much of an overlap (number of positions/nt) are required
our $REQUIRED_OVERLAP = 2;
# Length of ribosome shadow at beginning of UTR to mask (since miRNAs don't target right next to CDS)
our $BEG_UTR_MASK_LENGTH = 14;
# If $VERBOSE is non-zero, each Gene ID will be printed as it's processed. 
our $VERBOSE = 1;

our ($USAGE, $FILE_FORMATS, $GROUP_NUM, $MIR_FAM_ID, $LAST_UTR_ID);
our (@OUTPUT_THIS_GENE_THIS_MIR);
our (%MIR_ID_2_SEED, %MIR_ID_SPECIES, %MIR_TYPE_2_MATCH, %SPECIES_START_END, %SPECIES_START_END_2_MATCH, %SPECIES_TO_UTR, 
		%SPECIES_START_END_REMOVED, %SPECIES_START_END_2_MATCH_REMOVED, %GROUP_NUM_TO_SITE_TYPES, %GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST,
		%GROUP_TYPES_LIST_2_GROUP_TYPE, %SITE_TO_GROUP_NUM, %GROUP_NUM_TO_SPECIES, %GET_MATCH, %SITE_ID_2_SITE_TYPE, %SITE_ID_2_LENGTH,
		%SPECIES_START_END_MASKED);

# What types of seed matches should we look for? (1 ==> look for it; 0 ==> skip this one)
$GET_MATCH{1} = 1;	# 7mer-1a sites?
$GET_MATCH{2} = 1;	# 7mer-m8 sites?
$GET_MATCH{3} = 1;	# 8mer-1a sites?
$GET_MATCH{6} = 1;	# 6mer sites?
$GET_MATCH{4} = 0;	# 8mer-1u sites?
$GET_MATCH{5} = 0;	# 6mer-1a sites?

#######################  End of global variables  #######################

# Start program

getUsage();
getFileFormats();
our ($MIRNA_FILE, $UTR_FILE, $COORDS_FILE) = checkArguments();

# Get miRNA family data
readMiRNAs();

# Write conserved group to this file
open (COORDS, ">$COORDS_FILE") || die "Cannot open $COORDS_FILE for writing: $!";
# Print output file header
print COORDS "a_Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\tORF_overlap\n";

# Show how to summarize groups
getSiteTypeKeys();

# Go through UTR file, processsing one gene at a time
readUTRs();

############################  Subroutines  ############################


sub getUsage
{
	# Print command usage when arguments are missing 

	$USAGE = <<EODOCS;

	Description: Search for predicted miRNA targets
		     using the modified TargetScanS algorithm. 

	USAGE:
		$0 miRNA_file UTR_file PredictedTargetsOutputFile

	Required input files:
		miRNA_file    => miRNA families by species
		UTR_file      => Aligned UTRs		

	Output file:
		PredictedTargetsOutputFile    => Lists sites using alignment coordinates (MSA and UTR)

	For a description of input file formats, type
		$0 -h

	Author: George W. Bell, Bioinformatics and Research Computing
	Version: 7.0
	Copyright (c) The Whitehead Institute of Biomedical Research 

EODOCS
}


sub getFileFormats
{
	# Print input file formats with '-h' flag
	
	$FILE_FORMATS = <<EODOCS;

	** Required input files:
	
	1 - miRNA_file    => miRNA families by species
		
		contains three fields (tab-delimited):
			a. miRNA family ID/name
			b. seed region (7mer) for this miRNA
			c. semicolon-delimited list of species IDs in which this miRNA has been annotated
		ex:	   
		let-7/98	GAGGUAG	9606;10090;10116
		miR-127/127-3p	GAGGUAG	9606;10090
		
		Each miRNA family should be represented in a single line.
		
	2 - UTR_file      => Aligned UTRs		

		contains three fields (tab-delimited):
			a. Gene/UTR ID or name
			b. Species ID for this gene/UTR (must match ID in miRNA file)
			c. Aligned UTR or gene (with gaps from alignment)
		ex:
		BMP8B	9606	GUCCACCCGCCCGGC
		BMP8B	9615	-GUG--CUGCCCACC
		
		A gene will typically be represented on multiple adjacent lines.	

EODOCS
}


sub checkArguments
{
	# Check for input and output file arguments
	# Print info if there are any problems
	
	if ($ARGV[0] && $ARGV[0] eq "-h")
	{
		print STDERR "$USAGE";
		print STDERR "$FILE_FORMATS";
		exit (0);
	}
	elsif (! $ARGV[2])
	{
		print STDERR "$USAGE";
		exit (0);
	}
	elsif (! -e $ARGV[0])	# miRNA file not present
	{
		print STDERR "\nI can't find the file $ARGV[0]\n";
		print STDERR "which should contain the miRNA families by species.\n";
		exit;
	}
	elsif (! -e $ARGV[1])	# UTR file not present
	{
		print STDERR "\nI can't find the file $ARGV[1]\n";
		print STDERR "which should contain the aligned UTRs.\n";
		exit;
	}
	
	my $miRNAfile = $ARGV[0];
	my $UTRfile = $ARGV[1];
	my $coordsFile = $ARGV[2];
	
	if (-e $coordsFile)
	{
		print STDERR "Should I over-write $coordsFile [yes/no]? ";
		my $answer = <STDIN>;
		if ($answer !~ /^y/i)	{ exit; }
	}
	
	return ($miRNAfile, $UTRfile, $coordsFile);
}


sub readMiRNAs
{
	# Read user's file of miRNA families to use for target prediction

	my ($MIR_FAM_ID, $mirSeedRegion, $mirSpeciesIDlist);

	open (MIR_FAMILY_DATA, $MIRNA_FILE) || die "Cannot open $MIRNA_FILE for reading: $!";
	while (<MIR_FAMILY_DATA>)
	{
		# let-7/98	GAGGUAG	10090
	
		chomp;
		s/\r//g;	# For Windows and Mac
		
		###################  Public data format  ###################
		
		($MIR_FAM_ID, $mirSeedRegion, $mirSpeciesIDlist) = split (/\t/, $_);
		
		# Convert from RNA to DNA if needed
		$mirSeedRegion =~ s/T/U/gi;
		
		# Added by GB 28 Sep 2015
		$mirSeedRegion = uc($mirSeedRegion);
		$MIR_ID_2_SEED{$MIR_FAM_ID} = $mirSeedRegion;
		
		my @mirSpeciesIDlist = split /;/, $mirSpeciesIDlist;
		foreach my $mirSpeciesID (@mirSpeciesIDlist)
		{
			# Make sure we know which miRNA family is present in each species
			$MIR_ID_SPECIES{"${MIR_FAM_ID}::$mirSpeciesID"} = 1;
		}
	}

	# Get patterns for search for (one for each miRNA family)
	foreach $MIR_FAM_ID (sort keys %MIR_ID_2_SEED) 
	{
		get_seeds($MIR_ID_2_SEED{$MIR_FAM_ID}, $MIR_FAM_ID);
	}
}


sub get_seeds 
{
	# Get all types of seed matches for a given seed sequence

	my ($seedRegion, $MIR_FAM_ID) = @_;
	
	# Get the 7mer-m8 seed match (exactly the reverse complement)
	my $seed2 = $seedRegion;
	my $rseed2 = reverse($seed2);
	$rseed2 =~ tr/AUCG/UAGC/; 

	# Get the 6mer seed match (7mer-m8 minus the first position)
	my $rseed6 = $rseed2;
	$rseed6 =~ s/^.//;

	# Get the 6mer-1a seed match (6mer seed match minus the first position and add A at the end)
	my $rseed5 = $rseed6 . "A";
	$rseed5 =~ s/^.//;

	# Get the 7mer-1A seed match (6mer and add A at end)
	my $rseed1 = $rseed6 . "A";

	# Get the 8mer-A1 seed match (7mer-m8 seed match and add A at end)
	my $rseed3 = $rseed2 . 'A';

	# Get the 8mer-U1 seed match (7mer-m8 seed match and add U at end)
	my $rseed4 = $rseed2 . 'U';
	
	# print "SEED=$seedRegion\t7mer-1a=$rseed1\t7mer-m8=$rseed2\t8mer-1a=$rseed3\t8mer-1u=$rseed4\t6mer-1a=$rseed5\t6mer=$rseed6\n";
	
	# Make regex for searching by adding potential gaps between each pair of nts
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{1} = makeSeedMatchRegex($rseed1);
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{2} = makeSeedMatchRegex($rseed2);
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{3} = makeSeedMatchRegex($rseed3);
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{4} = makeSeedMatchRegex($rseed4);
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{5} = makeSeedMatchRegex($rseed5);
	$MIR_TYPE_2_MATCH{$MIR_FAM_ID}{6} = makeSeedMatchRegex($rseed6);	
}


sub makeSeedMatchRegex
{
	# Turn a seed match region into a Perl regular expression

	my $seedMatch = shift;
	my $seedMatchLength = length($seedMatch);
	my $seedMatchPattern = "";
	
	my @seedMatch = split '', $seedMatch;
	
	my $count = 0;	
	foreach my $seedMatchNt (@seedMatch) 
	{
		$count++;
		($count < $seedMatchLength) ? ($seedMatchPattern .= "$seedMatchNt-{0,}") : ($seedMatchPattern .= "$seedMatchNt");
	}
	return $seedMatchPattern;
}


sub getUTRcoords 
{
	###  Convert MSA coordinates to UTR coordinates
	
	my ($align, $end, $type) = ($_[0],$_[1],$_[2]);	
	
	my $utrBeg = substr($align, 0, $end);
	$utrBeg =~ s/\-//g;
	
	my $start = length($utrBeg); 
	$end = $start + $SITE_ID_2_LENGTH{$type} - 1;
	
	return ($start, $end); 
}


sub readUTRs
{
	# Initialize
	$LAST_UTR_ID = "";
	# This is the number (ID) of each group of sites; start with 1 and count up.
	$GROUP_NUM = 0;

	open (UTRS, $UTR_FILE) || die "Cannot open $UTR_FILE for reading: $!";
	while (<UTRS>)
	{
		# NM_031304       hg18    -       GGCCCCAC

		chomp;
		s/\r//g;	# For Windows and Mac  (When in doubt, convert input files to Unix format)

		if (! /^\s*$/)	# Ignore empty lines
		{
			# Public code format
			my ($utrID, $thisSpeciesID, $thisUTR) = split (/\t/, $_);

			if ($thisSpeciesID)	# Skip consensus sequence (if present)
			{
				# Convert from RNA to DNA if needed
				$thisUTR =~ s/T/U/gi;

				# Mask beginning of UTR since miRNA can't target UTR right next to CDS
				for (my $i = 0; $i < $BEG_UTR_MASK_LENGTH; $i++)
				{
					$thisUTR =~ s/[ACGTU]/N/i;
				}

				if ($utrID && $utrID ne $LAST_UTR_ID)
				{
					if ($VERBOSE)
					{
						# print STDERR "Processing $utrID\n";
						print "Processing $utrID\n";
					}

					# Look for sites in this gene (UTR) and process the results
					processUTRset();

					# Empty out these UTR-specific variables after finishing this set of UTRs
					%SPECIES_TO_UTR = ();
					%GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST = ();
					%GROUP_NUM_TO_SITE_TYPES = ();
					%GROUP_NUM_TO_SPECIES = ();
					%SITE_TO_GROUP_NUM = ();
				}

				# Add this UTR to the set
				$SPECIES_TO_UTR{$thisSpeciesID} = $thisUTR;

				$LAST_UTR_ID = $utrID; 			
			}
		}
	}
	# Get the last one
	processUTRset();
}


sub processUTRset
{
	# Look at each miRNA family
	foreach $MIR_FAM_ID (sort keys %MIR_ID_2_SEED) 
	{	
		# Do one species' UTR at a time
		foreach my $speciesIDthisUTR (sort keys %SPECIES_TO_UTR)
		{
			# Is this miRNA in this species?
			# If so [or if we want to look anyway], look for sites
			
			if ($FIND_SITES_ALL_SPECIES || $MIR_ID_SPECIES{"${MIR_FAM_ID}::$speciesIDthisUTR"})
			{
				# print STDERR "Looking for $MIR_FAM_ID sites in species ($speciesIDthisUTR)....\n";
			
				foreach my $matchType (keys %GET_MATCH)
				{
					if ($GET_MATCH{$matchType})
					{
						getMatches($MIR_FAM_ID, $speciesIDthisUTR, $matchType);
					}
				}

				###  Merge these types of sites when possible (new GB method)
				###  Start by dropping 6mer sites that are included by 7mer sites
				###  and then move to 7mer sites that are included in 8mer sites
				###  site types: # 6mer 6mer-1a 7mer-1a 7mer-m8 8mer-1a 8mer-1u
				
				findRemoveMatchSubsets();
			}
		}

		# If there are any hit(s) for this miRNA in any species, group each orthologous set of sites
		my $numSitesThisMirnaThisSpecies = 0;
		
		foreach my $site (sort keys %SPECIES_START_END)
		{
			$numSitesThisMirnaThisSpecies++;
		}

		if ($numSitesThisMirnaThisSpecies)
		{
			groupSitesThisGeneThisMiRNA();

			# Finish processing this gene/miR data 
			summarizePrintGroupsThisGeneThisMiRNA();
		}

		# Empty out these UTR+MIR-specific variables after finishing this $MIR_FAM_ID in this UTR
		@OUTPUT_THIS_GENE_THIS_MIR = ();
		%SPECIES_START_END = ();
		%SPECIES_START_END_2_MATCH = ();
		%SPECIES_START_END_REMOVED = ();
		%SPECIES_START_END_2_MATCH_REMOVED = ();
		%SPECIES_START_END_MASKED = ();
	}
}


sub getMatches
{
	# For a given seed sequence, return the positions and 
	# lengths of the matches in the sequence alignment (with gaps) 

	###  Corrected size of backtrack after match is found (21 Jul 2010)

	my ($MIR_FAM_ID, $speciesID, $matchType) = @_;
	my ($num, $start, $end, $matchedSubAlignment, $len);
	
	my $alignment = $SPECIES_TO_UTR{"$speciesID"};
	my $match = $MIR_TYPE_2_MATCH{$MIR_FAM_ID}{$matchType};
	
	# Convert to uppercase -- affects lower-case ORF-overlapping regions
	my $uppercaseAlignment = uc($alignment);
	
	while ($uppercaseAlignment =~ /$match/g) 
	{
		$end =  pos($uppercaseAlignment);
		$matchedSubAlignment = $&;
		$len = length($&);
		$start = $end - $len + 1;
		
		# Link site to site type
		$SPECIES_START_END{"${speciesID}::${start}::$end"} = $matchType;
		
		# Link site to actual match (sequence alignment region, possibly including gaps)
		$SPECIES_START_END_2_MATCH{"${speciesID}::${start}::$end"} = $matchedSubAlignment;
		
		# Compare original alignment to uppercase alignment; if they're different, flag as a masked region
		my $originalAlignedRegion = substr($alignment, $start, $len);
		my $ucAlignedRegion = substr($uppercaseAlignment, $start, $len);
		if ($ucAlignedRegion ne $originalAlignedRegion)
		{
			$SPECIES_START_END_MASKED{"${speciesID}::${start}::$end"} = 1;
			# print STDERR "${speciesID}::${start}::$end => $ucAlignedRegion != $originalAlignedRegion => masked!\n";
		}
		
		# pos($alignment) -= 5;	# but this may not be enough for matches with gaps
		# Backtrack far enough for matches with gaps (corrected 21 July 2010)
		pos($uppercaseAlignment) -= $len - 1;
	}
}


sub findRemoveMatchSubsets
{
	# Remove shorter matches that are a subset of longer matches
	# ex: 8mer-1a includes 7mer-m8 and 7mer-1a matches
	
	foreach my $siteThisUTRthisSpecies (keys %SPECIES_START_END)
	{
		my ($startPlusOne, $endMinusOne, $startPlusTwo);
		my ($species, $start, $end) = split (/::/, $siteThisUTRthisSpecies);
		my $gappedMatch = $SPECIES_START_END_2_MATCH{$siteThisUTRthisSpecies};
		
		if ($gappedMatch && $gappedMatch !~ /-/)
		{
			$startPlusOne = $start + 1;
			$startPlusTwo = $start + 2;
			$endMinusOne = $end - 1;
		}
		elsif ($gappedMatch)
		{
			my ($startPlusOneOffset, $startPlusTwoOffset, $endMinusOneOffset) = getSubsetCoords($gappedMatch);
			$startPlusOne = $start + $startPlusOneOffset;
			$startPlusTwo = $start + $startPlusTwoOffset;
			$endMinusOne = $end - $endMinusOneOffset;
			
			# print "Match with gaps: $gappedMatch ($siteThisUTRthisSpecies) [$startPlusOne, $startPlusTwo, $endMinusOne]\n";
		}
		
		my $dropSite = 0;

		if ($SPECIES_START_END{$siteThisUTRthisSpecies} && $SPECIES_START_END{$siteThisUTRthisSpecies} == 1)	# 7mer-1a
		{
			# Drop 6mer with same start position and 
			#      6mer-1a with same end position

			if ($GET_MATCH{6})	{ $dropSite = dropThisSite($species, $start, $endMinusOne); }	# 6mer
			if ($GET_MATCH{5})	{ $dropSite = dropThisSite($species, $startPlusOne, $end); }	# 6mer-1a			
		}
		elsif ($SPECIES_START_END{$siteThisUTRthisSpecies} && $SPECIES_START_END{$siteThisUTRthisSpecies} == 2)	# 7mer-m8
		{
			# Drop 6mer with same end position
			
			if ($GET_MATCH{6})	{ $dropSite = dropThisSite($species, $startPlusOne, $end); }	# 6mer
		}
		elsif ($SPECIES_START_END{$siteThisUTRthisSpecies} && $SPECIES_START_END{$siteThisUTRthisSpecies} == 3)	# 8mer-1a
		{
			# Drop 7mer-m8 with same starting position and
			#      7mer-1a with same ending position and 
			#      6mer-1a with same ending position and
			#      6mer starting one position later

			if ($GET_MATCH{2})	{ $dropSite = dropThisSite($species, $start, $endMinusOne); }	# 7mer-m8
			if ($GET_MATCH{1})	{ $dropSite = dropThisSite($species, $startPlusOne, $end); }	# 7mer-1a
			if ($GET_MATCH{5})	{ $dropSite = dropThisSite($species, $startPlusTwo, $end); }	# 6mer-1a
			if ($GET_MATCH{6})	{ $dropSite = dropThisSite($species, $startPlusOne, $endMinusOne); }	# 6mer			
		}
		elsif($SPECIES_START_END{$siteThisUTRthisSpecies} && $SPECIES_START_END{$siteThisUTRthisSpecies} == 4)	# 8mer-1u
		{
			# Drop 7mer-m8 with same starting position and 
			#      6mer starting one position later

			if ($GET_MATCH{2})	{ $dropSite = dropThisSite($species, $start, $endMinusOne); }	# 7mer-m8
			if ($GET_MATCH{6})	{ $dropSite = dropThisSite($species, $startPlusOne, $endMinusOne); }	# 6mer
		}
	}
}


sub getSubsetCoords
{
	# Given an alignment with gaps, we need to identify the positions of the
	# second, third, and second-to-last nucleotides (so we can drop the appropriate match subsets)
	
	my $alignment = shift;
	my ($matchPos, $startPlusOneOffset, $startPlusTwoOffset, $endMinusOneOffset);
	
	$alignment =~ /^[^-]-*([^-])/g;
	$matchPos = pos($alignment);
	$startPlusOneOffset = $matchPos - 1;
	# print "startPlusOne match in $alignment at $matchPos ($1), so offset is $startPlusOneOffset\n";
	pos($alignment) = 0;

	$alignment =~ /^[^-]-*[^-]-*([^-])/g;
	$matchPos = pos($alignment);
	$startPlusTwoOffset = $matchPos - 1;
	# print "startPlusTwo match in $alignment at $matchPos ($1), so offset is $startPlusTwoOffset\n";
	pos($alignment) = 0;

	$alignment =~ /([^-])-*[^-]$/g;
	$endMinusOneOffset = length($&) - 1;
	# print "endMinusOne match in $alignment ($1), so offset is $endMinusOneOffset\n\n";

	return ($startPlusOneOffset, $startPlusTwoOffset, $endMinusOneOffset);
}


sub dropThisSite
{
	# Drop a site that is a subset of another site
	my ($species, $start, $end) = @_;
	
	if ($SPECIES_START_END{"${species}::${start}::$end"})
	{
		# Keep record of what we deleted
		$SPECIES_START_END_REMOVED{"${species}::${start}::$end"} = $SPECIES_START_END{"${species}::${start}::$end"};  
		$SPECIES_START_END_2_MATCH_REMOVED{"${species}::${start}::$end"} = $SPECIES_START_END_2_MATCH{"${species}::${start}::$end"};

		delete $SPECIES_START_END{"${species}::${start}::$end"};  
		delete $SPECIES_START_END_2_MATCH{"${species}::${start}::$end"};
		
		# print "We're deleting ${species}::${start}::$end\n";

		return 1;
	}
	elsif ($SPECIES_START_END_REMOVED{"${species}::${start}::$end"})
	{
		# We already deleted this site
		# print "We already deleted ${species}::${start}::$end\n";
		return 1;
	}
	else	# This site doesn't exist
	{
		# Was it already deleted, or is there a gap that requires adjustment???	
		return 0;
	}
}


sub groupSitesThisGeneThisMiRNA
{
	# Group miRNA sites that overlap in alignment
	
	# Initialize
	%SITE_TO_GROUP_NUM = ();
	my $numOverlapNt;
	
	# Can we drop this variable???  It doesn't seem to be used.
	my %pairToDistance;

	############### Check for position overlap between sites

	# Do an all vs. all comparison to identify overlaps
	foreach my $site1 (sort keys %SPECIES_START_END)
	{
		my ($site1Species, $site1Start, $site1End) = split (/::/, $site1);
	
		foreach my $site2 (sort keys %SPECIES_START_END)
		{
			my ($site2Species, $site2Start, $site2End) = split (/::/, $site2);
			
			# Skip comparison of same-species sites
			if ($site1Species ne $site2Species)
			{
				############  Choose combinations to give overlap  ############

				# Same start and end
				if ($site1Start eq $site2Start && $site1End eq $site2End)
				{
					groupThisPair($site1, $site2);
					
					$pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;
				}
				# Same start
				elsif ($site1Start eq $site2Start)
				{
					groupThisPair($site1, $site2);
					
					if ($site1End > $site2End)
					{ $pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;}
					else { $pairToDistance{"$site1 $site2"} = $site2End - $site2Start + 1;}
				}			
				# Same end
				elsif ($site1End eq $site2End)
				{
					groupThisPair($site1, $site2);
					
					if ($site1Start < $site2Start)
					{ $pairToDistance{"$site1 $site2"} = $site1End - $site1Start + 1;}
					else { $pairToDistance{"$site1 $site2"} = $site2End - $site2Start + 1;}
					
				}
				# Offset one direction
				#     xxxxxxx
				#    xxxxxxx
				elsif ($site1Start > $site2Start && $site1Start <= $site2End)
				{
					$numOverlapNt = $site2End - $site1Start + 1;
					if ($numOverlapNt >= $REQUIRED_OVERLAP)
					{
						groupThisPair($site1, $site2);
						$pairToDistance{"$site1 $site2"} = $numOverlapNt;
					}				
				}
				# Offset other direction
				#    xxxxxxx
				#     xxxxxxx
				elsif ($site1End >= $site2Start && $site1End < $site2End)
				{
					$numOverlapNt = $site1End - $site2Start + 1;
					if ($numOverlapNt >= $REQUIRED_OVERLAP)
					{
						groupThisPair($site1, $site2);
						$pairToDistance{"$site1 $site2"} = $numOverlapNt;
					}
				}
				# One within the other (with gaps)
				#      xxxxxxx          xxxxxxxxx
				#     xxxxxxxxx          xxxxxxx
				elsif ( ($site1Start > $site2Start && $site1End < $site2End) ||
					($site2Start > $site1Start && $site2End < $site1End) )
				{
					groupThisPair($site1, $site2);
					$pairToDistance{"$site1 $site2"} = $numOverlapNt;
				}				
			}
		}
	}
	
	foreach my $thisSite (sort keys %SPECIES_START_END)
	{
		my $annotated;
		
		my @siteAllInfo = split (/::/, $thisSite);
		
		my $speciesThisSite = $siteAllInfo[0];
		
		# print "Site is *$thisSite* and has GROUP_NUM $SITE_TO_GROUP_NUM{$thisSite}\n";
	
		# This site is a group of 1, so no group info yet
		if (! $SITE_TO_GROUP_NUM{$thisSite})
		{
			# If this group hasn't yet been assigned a number, give it one.
			$GROUP_NUM++;
			# push @groupNumThisGeneThisMir, $GROUP_NUM;
			$SITE_TO_GROUP_NUM{$thisSite} = $GROUP_NUM;
			
			# print "Now assigned $thisSite to GROUP_NUM $GROUP_NUM\n";
		}
		
		if (! $GROUP_NUM_TO_SITE_TYPES{$SITE_TO_GROUP_NUM{$thisSite}})
		{
			# Start a list of site types for this group
			$GROUP_NUM_TO_SITE_TYPES{$SITE_TO_GROUP_NUM{$thisSite}} = "$SPECIES_START_END{$thisSite}";
		}
		else
		{
			# Add to the list of site types for this group
			$GROUP_NUM_TO_SITE_TYPES{$SITE_TO_GROUP_NUM{$thisSite}} .= ";$SPECIES_START_END{$thisSite}";
		}
		
		# Make a list of species in which a site type is found (in this group)  
		$GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{$SPECIES_START_END{$thisSite}} .= "$speciesThisSite ";
						
		###  If a wide site is present, its subset sites are also present
		if    ($SPECIES_START_END{$thisSite} == 1)	# 7mer-1a
		{
			if ($GET_MATCH{6})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{6} .= "$speciesThisSite "; }	# 6mer
			if ($GET_MATCH{5})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{5} .= "$speciesThisSite "; }	# 6mer-1a
		}
		elsif ($SPECIES_START_END{$thisSite} == 2)	# 7mer-m8
		{
			if ($GET_MATCH{6})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{6} .= "$speciesThisSite "; }	# 6mer
		}
		elsif ($SPECIES_START_END{$thisSite} == 3)	# 8mer-1a
		{
			if ($GET_MATCH{1})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{1} .= "$speciesThisSite "; }	# 7mer-1a
			if ($GET_MATCH{2})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{2} .= "$speciesThisSite "; }	# 7mer-m8
			if ($GET_MATCH{5})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{5} .= "$speciesThisSite "; }	# 6mer-1a
			if ($GET_MATCH{6})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{6} .= "$speciesThisSite "; }	# 6mer
		}		
		elsif ($SPECIES_START_END{$thisSite} == 4)	# 8mer-1u
		{
			if ($GET_MATCH{2})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{2} .= "$speciesThisSite "; }	# 7mer-m8
			if ($GET_MATCH{6})	{ $GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$SITE_TO_GROUP_NUM{$thisSite}}{6} .= "$speciesThisSite "; }	# 6mer
		}
		###

		# Given the MSA coords, get the corresponding UTR coords
		# Correct 8mer OBOB: Nov 26, 2007
		my ($utrStart, $utrEnd) = getUTRcoords($SPECIES_TO_UTR{$siteAllInfo[0]}, $siteAllInfo[1], $SPECIES_START_END{$thisSite});
	
		# Link each group to the species within it
		$GROUP_NUM_TO_SPECIES{$SITE_TO_GROUP_NUM{$thisSite}} .= "$siteAllInfo[0];";
	
		# Is this miRNA annotated in this species? 
		if (! $MIR_ID_SPECIES{"${MIR_FAM_ID}::$speciesThisSite"})
		{
			$annotated = " ";
		}
		else
		{
			$annotated = "x";
		}
		
		push @OUTPUT_THIS_GENE_THIS_MIR, "$LAST_UTR_ID\t$MIR_FAM_ID\t$speciesThisSite\t$siteAllInfo[1]\t$siteAllInfo[2]\t$utrStart\t$utrEnd\t$SITE_TO_GROUP_NUM{$thisSite}\t$SPECIES_START_END{$thisSite}\t$annotated";
	}
}


sub groupThisPair
{
	# Take a pair of overlapping sites (same or different species)
	# If one of them has been assigned to a $GROUP_NUM, assign the other one
	# If neither has been assigned yet, increment $GROUP_NUM and assign them both

	my ($site1, $site2) = @_;

	if (! $SITE_TO_GROUP_NUM{$site1} || ! $SITE_TO_GROUP_NUM{$site2})
	{
		if ($SITE_TO_GROUP_NUM{$site1})	# Site 1 already part of a group
		{
			# Set $site2 to the same group as $site1
			$SITE_TO_GROUP_NUM{$site2} = $SITE_TO_GROUP_NUM{$site1};
		}
		elsif ($SITE_TO_GROUP_NUM{$site2})	# Site 2 already part of a group
		{
			# Set $site1 to the same group as $site2
			$SITE_TO_GROUP_NUM{$site1} = $SITE_TO_GROUP_NUM{$site2};
		}
		else
		{
			# Increment the group number
			$GROUP_NUM++;
			$SITE_TO_GROUP_NUM{$site1} = $GROUP_NUM;
			$SITE_TO_GROUP_NUM{$site2} = $GROUP_NUM;
		}
	}
}


sub summarizePrintGroupsThisGeneThisMiRNA
{
	# print "Summarizing & printing sites for $LAST_UTR_ID ($MIR_FAM_ID)...\n";

	my (@species,@uniqueSpecies, @speciesThisGroup);
	my (%groupNumToSiteTypesList, %groupToInfo);
	my ($uniqueSiteTypesList, $uniqueSiteTypesNamesList, $isMasked);
	%groupToInfo = ();

	foreach my $groupNum (sort {$a <=> $b} keys %GROUP_NUM_TO_SITE_TYPES ) 
	{
		@speciesThisGroup = makeListNonRedundant($GROUP_NUM_TO_SPECIES{$groupNum}, ";");
	
		if ($GROUP_NUM_TO_SITE_TYPES{$groupNum})
		{
			my @siteTypes = split (/;/, $GROUP_NUM_TO_SITE_TYPES{$groupNum});
			
			# Make this site type list unique (shouldn't be necessary) and sorted 5.2.07
			my @uniqueSiteTypes = makeListNonRedundant($GROUP_NUM_TO_SITE_TYPES{$groupNum}, ";");
			
			# Convert site types list into site names
			$uniqueSiteTypesNamesList = "";
			foreach my $siteType (@uniqueSiteTypes)
			{
				if (! $uniqueSiteTypesNamesList)
				{
					$uniqueSiteTypesNamesList = $SITE_ID_2_SITE_TYPE{$siteType};
				}
				else
				{
					$uniqueSiteTypesNamesList .= "+$SITE_ID_2_SITE_TYPE{$siteType}";
				}
			}

		}
		else	# There's a problem here
		{
			$uniqueSiteTypesList = "";
		}
		
		$groupNumToSiteTypesList{$groupNum} = $uniqueSiteTypesNamesList;
		
		$groupToInfo{$groupNum} = "@speciesThisGroup";		
	}

	###  Sort array by "group number" field
	@OUTPUT_THIS_GENE_THIS_MIR = map  { $_->[0] }
	sort {
			$a->[1] <=> $b->[1]  # first of the selected fields
		 }
	map  { [ $_, (split /\t/)[7] ] } # select desired fields of a tab-delimited line in array
	@OUTPUT_THIS_GENE_THIS_MIR; 
	
	foreach my $dataOneSiteThisGeneThisMir (@OUTPUT_THIS_GENE_THIS_MIR)
	{
		my @f = split (/\t/, $dataOneSiteThisGeneThisMir);
		my $groupNumThisSite = $f[7];
		my $siteTypeThisSite = $f[8];
		my $groupType;
		# Replace site type ID with site type name
		$f[8] = $SITE_ID_2_SITE_TYPE{$siteTypeThisSite};
		$dataOneSiteThisGeneThisMir = join "\t", @f;
		
		if ($GROUP_TYPES_LIST_2_GROUP_TYPE{$groupNumToSiteTypesList{$groupNumThisSite}})
		{
			$groupType = $GROUP_TYPES_LIST_2_GROUP_TYPE{$groupNumToSiteTypesList{$groupNumThisSite}};
		}
		else	# Unexpected set of groups
		{
			$groupType = $groupNumToSiteTypesList{$groupNumThisSite};
		}
		
		$dataOneSiteThisGeneThisMir .= "\t$groupType";
		$dataOneSiteThisGeneThisMir .= "\t$groupToInfo{$groupNumThisSite}";
		
		# Add the species list for this site type (but only if site type ne group type)
		if ($groupType ne $SITE_ID_2_SITE_TYPE{$siteTypeThisSite})
		{
			my @nrSpeciesThisSubset = makeListNonRedundant($GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST{$groupNumThisSite}{$siteTypeThisSite}, " ");
			my $nrSpeciesThisSubset = join " ", @nrSpeciesThisSubset;
			$dataOneSiteThisGeneThisMir .= "\t$nrSpeciesThisSubset";
		}
		else
		{
			$dataOneSiteThisGeneThisMir .= "\t";
		}
		
		$isMasked = $SPECIES_START_END_MASKED{"$f[2]::$f[3]::$f[4]"} ? 1 : 0;
		$dataOneSiteThisGeneThisMir .= "\t$isMasked";	

		print COORDS "$dataOneSiteThisGeneThisMir\n";
	}
}


sub makeListNonRedundant
{
	# Convert a list string into a non-redundant array
	my($list, $sep) = @_;

	my @values = split ($sep, $list);

	my %uniqueList = ();
	%uniqueList = map { $_ => 1 } @values;
	my @uniqueList = keys %uniqueList;		
	@uniqueList = sort @uniqueList;
	
	# If entries are numbers
	@uniqueList = sort {$a <=> $b} @uniqueList;
	# If entries are not numbers
	# @uniqueList = sort @uniqueList;

	return @uniqueList;
}


sub getSiteTypeKeys
{
	# Convert site type ID into name
	$SITE_ID_2_SITE_TYPE{1} = "7mer-1a";
	$SITE_ID_2_SITE_TYPE{2} = "7mer-m8";
	$SITE_ID_2_SITE_TYPE{3} = "8mer-1a";
	$SITE_ID_2_SITE_TYPE{4} = "8mer-1u";
	$SITE_ID_2_SITE_TYPE{5} = "6mer-1a";
	$SITE_ID_2_SITE_TYPE{6} = "6mer";
	
	$SITE_ID_2_LENGTH{1} = 7;
	$SITE_ID_2_LENGTH{2} = 7;
	$SITE_ID_2_LENGTH{3} = 8;
	$SITE_ID_2_LENGTH{4} = 8;
	$SITE_ID_2_LENGTH{5} = 6;
	$SITE_ID_2_LENGTH{6} = 6;
}
