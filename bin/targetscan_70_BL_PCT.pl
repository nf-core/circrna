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
# Version: 6.0 (31 October 2011)
# Version: 7.0 (30 September 2014 - February 2015)
#
# Comment: This program calculates branch length and probability of conserved targeting
#          for miRNA target sites as predicted by TargetScan 6
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70
#
# Version 6 is identical to version 5 except the three site type names were modified.
#
# Version 7: Branch length is calculated only for sites of miRNAs present in input file miRNA_file
#            expand the set of species
#            update parameters           
#            process 6mers (if in input file)
#            use masked field (if present) to set PCT to "NA"
#
#######################################################################

# Calculate BL and then PCT for each site in a group,
# processing one group (several lines of data) 
# Input file is output from targetscan_70_BaRC.pl
# with new format including "Species_in_this_group_with_this_site_type" column

# Input file must first be sorted by Group_num (column 8), like
# sort -k8,8n UTRs.00.23species.targets.txt > UTRs.00.23species.targets.sort.txt

# If site type == group type, use group type and group species for BL calculation
# If site type != group type, try all possible group types for this species,
# using appropriate set of species in which this type is conserved. 
# As a result, each row will have 0 (for single-species sites) BLs to 3 BLs
# (where the reference species is an 8mer but the group also contains 7mer-m8 and 7mer-1a sites.
# Select the BL and group type to get the best group with a BL that meets the conservation threshold. 
#
#

# Need this to read species trees, so BioPerl needs to be installed
use Bio::TreeIO;

# e -- used to calculation PCT
$eConstant = 2.71828182845904523536;

# PCT data files with S:B data for each highly conserved seed for a series of branch lengths
# these have been corrected so there are no negative BG values
$pathToTreesParams = "PCT_parameters";	# This dir also contains 10 trees, one for each bin
$PCT_data_8mer    = "$pathToTreesParams/8mer_PCT_parameters.txt";
$PCT_data_7mer_m8 = "$pathToTreesParams/7mer_m8_PCT_parameters.txt";
$PCT_data_7mer_a1 = "$pathToTreesParams/7mer_1a_PCT_parameters.txt";

# Names of site types in predicted targets file
$siteType3name = "8mer-1a";
$siteType2name = "7mer-m8";
$siteType1name = "7mer-1a";
$siteType4name = "6mer";

# If the BL is >= these thresholds, it's "conserved" (for TS7 trees)
# For TargetScanHuman
$siteTypeToConsThreshold{$siteType3name} = 1.8;
$siteTypeToConsThreshold{$siteType2name} = 2.8;
$siteTypeToConsThreshold{$siteType1name} = 3.6;
$siteTypeToConsThreshold{$siteType4name} = 100; # all should be nonconserved
# For TargetScanMouse
# $siteTypeToConsThreshold{$siteType3name} = 0.6;
# $siteTypeToConsThreshold{$siteType2name} = 1.8;
# $siteTypeToConsThreshold{$siteType1name} = 2.5;
# $siteTypeToConsThreshold{$siteType4name} = 100; # all should be nonconserved

# What are all the nodes in the tree?
# These are NCBI Taxonomy IDs included in the analysis
@all_orgs = qw ( 9606 9598 9595 9601 61853 9544 9541 9557 60711 9483 27679 30611 246437 43179 51337 79684 10029 10036 10090 10116 10181 10141 34839 10160 9986 9978 9823 30538 419612 9739 9733 59538 9913 9940 9925 9796 9807 9685 9615 9669 9646 9708 9713 9402 132908 225400 59463 29078 9365 42254 143302 9785 28737 127582 185453 9371 1230840 9361 13616 9305 9315 9258 345164 8954 59894 44394 48883 59729 181119 13146 241585 176014 8932 8839 9031 9103 8496 8469 8478 13735 55534 28377 8364 7897 );

$NULL_VALUE = "NA";

# Note: We can safely ignore the several errors from Bio/Tree/Node.pm at the end of the analysis.

###########################################  End of main constants  ###########################################

# Check for input arguments
if (! $ARGV[2])
{
	print STDERR "\nCalculate branch length and PCT for miRNA targets predicted by TargetScan\n";
	print STDERR "using a BL tree specific to each gene (UTR) bin\n";
	print STDERR "\nUSAGE: $0 miRNA_file predicted_targets UTR_bin_info > out\n";
	print STDERR "Ex: $0 miR_Family_info_sample.txt targetscan_70_output.txt UTRs_median_BLs_bins.txt > targetscan_70_output.BL_PCT.txt\n\n";
	print STDERR "Branch length is calculated only for sites of miRNAs present in input file miRNA_file.\n";
	print STDERR "PCT is only relevant for highly conserved miRNA families.\n\n";
	exit;
}

# TargetScan code input file
$miRNAfile = $ARGV[0];
# TargetScan code output file
$predictedTargets = $ARGV[1];
# Output from get_median_BLs_bins_from_alignments.pl
$geneBinFile = $ARGV[2];

# Has the $predictedTargets file been previously sorted?  If so, no need to re-sort.
# If this file includes the word "sort", we're assuming it's already been sorted.
if ($predictedTargets !~ /sort/)
{
	# Targets need to be sorted by group number
	# This requires a Linux-style 'sort' command on the system
	
	$predictedTargets = sortPredictedTargets($predictedTargets);
}
else
{
	print STDERR "\nWe're assuming $predictedTargets has already been sorted, so we're not sorting it.\n";
	print STDERR "If it hasn't already been sorted, please remove \"sort\" from the file name.\n\n";
}

# Get miRNA family ==> seed conversion
readMiRNAs();

# Get bin for each gene (UTR)
readBLbins();

# Read all files (and do a little math) needed to get PCT values in memory
readAllPCTdata();

# Go through the predicted sites, getting a list of species for this site 
# and then, using a tree specific to the bin, calculate the branch length

open (TARGET_SITES_IN, $predictedTargets) || die "Cannot open $predictedTargets for reading: $!";
while (<TARGET_SITES_IN>)
{
	# Read the predicted targets file.
	
	# a_Gene_ID	miRNA_family_ID	species_ID	MSA_start	MSA_end	UTR_start	UTR_end	Group_num	Site_type	miRNA in this species	Group_type	Species_in_this_group	Species_in_this_group_with_this_site_type
	# CDC2L6	miR-1/206	9606	2827	2833	2588	2594	1	7mer-1a	x	7mer-1a	9606
	# CDC2L6	miR-1/206	9615	2177	2183	2053	2059	2	7mer-1a	x	7mer-1a	9615
	
	chomp;
	my $targetLine = $_;
	
	if ($. > 1)	# Skip header line
	{
		my @f = split (/\t/, $targetLine);

		$miRNA_family = $f[1];
		$groupNum = $f[7];
		$siteType = $f[8];
		$groupSpeciesList = $f[11];
		$groupSiteTypeSpeciesList = $f[12];

		# Starting on a new group
		if ($prevGroupNum && $groupNum ne $prevGroupNum)
		{
			# Analyze the previous group
			processPreviousGroup();		
		}

		# Skip miRNA families not in list from miRNA_file
		if ($groupNum =~ /\d+/ && $getBLSthisFamily{$miRNA_family})	
		{
			# Make a hash of all lines in this group 
			$lineNum2data{$.} = $targetLine;

			# Check that we have the correct group type
			push @sitesTypesThisGroup, $siteType;

			if ($groupSiteTypeSpeciesList)
			{
				$siteType2speciesList{$siteType} = $groupSiteTypeSpeciesList;
			}
			else
			{
				$siteType2speciesList{$siteType} = $groupSpeciesList;
			}
		}
		else
		{
			# Print input row, adding NA for BLS and PCT fields
			$f[11] = "NA";
			$f[12] = "NA";
			print join "\t", @f, "\n";
		}
		$prevGroupNum = $groupNum;
	}
	else
	{
		print "Gene_ID\tmiRNA_family_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tBranch length score\tPct\tConserved\n";
	}
}

# Get the last one
processPreviousGroup();

##########################################################  Subroutines  ##########################################################

sub sortPredictedTargets
{
	# This sorting requires a Linux-style sort on our computer
	# If it's not available, sort using another method and names the sorted file by ending your file
	# by foo.sort.txt instead of foo.txt
	# Use foo.txt as the input file, but we'll see if the sorted version exists and use if 
	# instead of trying to sort

	$predictedTargets = $_[0];
	
	$predictedTargetsSorted = $predictedTargets;
	
	# Make new name for sorted file
	if ($predictedTargets =~ /\.txt$/)
	{
		$predictedTargetsSorted =~ s/\.txt$/.sort.txt/;
	}
	else
	{
		$predictedTargetsSorted .= ".sort.txt";
	}
	
	# If sorted version doesn't exist, create it
	if (! -e $predictedTargetsSorted)
	{
		$CMD = "sort -k8,8n $predictedTargets > $predictedTargetsSorted";
	
		print STDERR "Sorting $predictedTargets by group number and creating $predictedTargetsSorted\n";
	
		`$CMD`;
		
		print STDERR "Sorting [ $CMD ] is done\n";
	}
	else	# If sorted version does exist, use it instead.
	{
		print STDERR "Using $predictedTargetsSorted as input file.  If this is wrong, delete this file and try again.\n";
	}
	
	return $predictedTargetsSorted;
}

sub readMiRNAs
{
	# Read the miRNA file to link a miRNA name to its seed region.
	# The third field of the file is ignored.
	# New for TS7: Include only miRNA families that are conserved or broadly conserved (so NOT poorly conserved ones)

	my ($mirFamID, $mirSeedRegion);

	open (MIR_FAMILY_DATA, $miRNAfile) || die "Cannot open $miRNAfile for reading: $!";
	while (<MIR_FAMILY_DATA>)
	{
		# let-7/98	GAGGUAG	10090
	
		chomp;
		s/\r//g;	# For Windows and Mac

		my @f = split (/\t/, $_);
		
		$mirFamID = $f[0];
		$mirSeedRegion = $f[1];
		
		# Convert from RNA to DNA if needed
		$mirSeedRegion =~ s/T/U/gi;

		$mirID2seed{$mirFamID} = $mirSeedRegion;
		
		# New for TS7: keep track of families for which we want to calculate BLS
		$getBLSthisFamily{$mirFamID} = 1;
		$getBLSthisFamily{$mirSeedRegion} = 1;
	}
}

sub readBLbins
{
	# Open a file (generated by targetscan_60_BL_bins.pl)
	# with gene IDs in first column and bin numbers (1 - 10) in third column
	# The second column is ignored.

	open (GENE2BIN, $geneBinFile) || die "Cannot open $geneBinFile for reading: $!";
	while (<GENE2BIN>)
	{
		# CDC2L6	0.63984	1
		# FNDC3A	1.15171	2

		chomp;

		my @f = split (/\t/, $_);

		# Simplify input file (18 Sep 2014)
		$refseq2bin{$f[0]} = $f[2]; 
	}
	close (GENE2BIN);
}

sub readAllPCTdata
{
	# Read background data (3 files) needed to calculate PCT for highly conserved miRNA families
	# These were generated by Robin Friedman, MIT.

	push @PCT_data_files, $PCT_data_8mer;
	push @PCT_data_files, $PCT_data_7mer_m8;
	push @PCT_data_files, $PCT_data_7mer_a1;

	# Link data files to site types
	$sbFile2SiteType{$PCT_data_8mer} = $siteType3name;
	$sbFile2SiteType{$PCT_data_7mer_m8} = $siteType2name;
	$sbFile2SiteType{$PCT_data_7mer_a1} = $siteType1name;

	foreach my $sbFile (@PCT_data_files)
	{
		open (SB_DATA, $sbFile) || die "Cannot open $sbFile for reading: $!\n";

		# print STDERR "Reading PCT data from $sbFile\n";

		my $siteType = $sbFile2SiteType{$sbFile};

		my $lineNum = 0;

		while(<SB_DATA>)
		{
			chomp;
			my @f = split(/\t/, $_);

			if ($lineNum)	# Data occurs after first (header) line
			{
				my $miRNAfamily = $f[0];
				
				$thisMiRNAhasPCT{$miRNAfamily} = 1;

				# Get coefficients b0, b1, b2, and b3
				for (my $i = 1; $i <= 4; $i++)	# Parameters in columns 2 - 5
				{
					$familyPlusType2coeff{$miRNAfamily}{$siteType}{$i - 1} = $f[$i];	# Parameters in columns 2 - 5
				}
			}
			$lineNum++;
		}
		close(SB_DATA);
	}
}

sub processPreviousGroup
{
	# Process one group of miRNA sites (overlapping sites of one or more species) at a time

	# Check and update the group assignment  
	my $recheckedGroupType = getGroupTypeAnotherWay(@sitesTypesThisGroup);

	foreach $lineNum (sort {$a <=> $b} keys %lineNum2data)
	{
		$groupLine = $lineNum2data{$lineNum};
	
		chomp($groupLine);
		
		my @f = split (/\t/, $groupLine);

		my $refseq = $f[0];
		
		my $miRNAfamily;
		
		if ($mirID2seed{$f[1]})
		{
			$miRNAfamily = $mirID2seed{$f[1]};	# This conversion comes from $miRNAfile
		}
		else
		{
			print STDERR "No seed for $f[1]\n";
			$miRNAfamily = "UNKNOWN";
		}
		my $groupNum = $f[7];
		my $siteType = $f[8];
		my $groupType = $f[10];
		my $groupSiteTypeSpeciesList = $f[12];
		my $blAll = "";
		my $pctAll = "";
		my $pctMax = "";
		my @pct = ();
		my $pct = "";
		my $isMasked = $f[13];
		
		if ($groupType ne $recheckedGroupType)
		{
			$groupType = $recheckedGroupType;
			$f[10] = $recheckedGroupType;
		}
		
		$bl_8mer = $NULL_VALUE;
		$bl_m8 = $NULL_VALUE;
		$bl_1a = $NULL_VALUE;
		
		if ($siteType eq $groupType || $siteType ne "8mer")
		{
			# Possible combinations:
			# Site_type       Group_type
			# 8mer    8mer
			# m8      m8		
			# 1a      1a		
			# 1a      8mer+1a
			# 1a      8mer+m8+1a
			# 1a      m8+1a
			# m8      8mer+m8
			# m8      8mer+m8+1a
			# m8      m8+1a
		
			if ($siteType2speciesList{$siteType})
			{
				$blSelected = lookAtBranchLength($refseq, $siteType2speciesList{$siteType});
				# $blSelected = ($siteType ne "6mer")? lookAtBranchLength($refseq, $siteType2speciesList{$siteType}) : 0;
			}
			else
			{
				print STDERR "No list for $groupNum site type of $_\n";
			}
			
			$groupTypeForThisSite = $siteType;

			# Use BL and site type to determine whether this site is conserved
			# This may be different for different sites in the group

			$consScore = selectConservationForThisSite($blSelected, $groupTypeForThisSite);
		}

		elsif ($siteType eq "8mer-1a" && $siteType2speciesList{"8mer-1a"})
		{
			# Possible combinations:
			# Site_type       Group_type
			# 8mer    8mer+1a
			# 8mer    8mer+m8
			# 8mer    8mer+m8+1a
		
			if ($groupType =~ /8mer-1a/ && $siteType2speciesList{"8mer-1a"})
			{
				$bl_8mer = lookAtBranchLength($refseq, $siteType2speciesList{"8mer-1a"});
			}
			if ($groupType =~ /7mer-m8/ && $siteType2speciesList{"7mer-m8"})
			{
				# Possible combinations:
				# Site_type       Group_type
				# 8mer    8mer+m8
				# 8mer    8mer+m8+1a
	
				$bl_m8 = lookAtBranchLength($refseq, $siteType2speciesList{"m8"});
			}
			if ($groupType =~ /7mer-1a/ && $siteType2speciesList{"7mer-1a"})
			{
				# Possible combinations:
				# Site_type       Group_type
				# 8mer    8mer+1a
				# 8mer    8mer+m8+1a
			
				$bl_1a = lookAtBranchLength($refseq, $siteType2speciesList{"1a"});
			}
			
			$blAll = "${bl_8mer};${bl_m8};${bl_1a}";
			
			$groupTypePlusBL = selectGroupTypeBLforThisSite($blAll);
			
			($groupTypeForThisSite, $blSelected) = split / /, $groupTypePlusBL;
			
			$consScore = selectConservationForThisSite($blSelected, $groupTypeForThisSite);
		}
		# Check here is this site is in a masked UTR region (overlapping an ORF).  If so, set PCT to $NULL_VALUE.
		if ($thisMiRNAhasPCT{$miRNAfamily} && ! $isMasked)
		{
			if (! $blAll)	# Only one possible BL for this site
			{
				$pct = getPCT($groupTypeForThisSite, $miRNAfamily, $blSelected);
				# $siteTypeOfPctMax = $groupTypeForThisSite;
			}
			else	# Want to try getting three PCTs and choose highest one
			{
				@blAll = split /;/, $blAll;
				
				# Initialize
				$pct[0] = $pct[1] = $pct[2] = $NULL_VALUE;
				
				if ($blAll[0] ne $NULL_VALUE && $blAll[0] > 0) # 8mer site
				{
					$pct[0] = getPCT($siteType3name, $miRNAfamily, $blAll[0]);
					$pctMax = $pct[0];
					$siteTypeOfPctMax = $siteType3name;
				}
				if ($blAll[1] ne $NULL_VALUE && $blAll[1] > 0) # 7mer-m8 site
				{  
					$pct[1] = getPCT($siteType2name, $miRNAfamily, $blAll[1]);
					if (! $pctMax || $pct[1] > $pctMax)
					{
						$pctMax = $pct[1];
						$siteTypeOfPctMax = $siteType2name;
					}
				}	
				if ($blAll[2] ne $NULL_VALUE && $blAll[2] > 0) # 7mer-1a site
				{  
					$pct[2] = getPCT($siteType1name, $miRNAfamily, $blAll[2]);
					if (! $pctMax || $pct[2] > $pctMax)
					{
						$pctMax = $pct[2];
						$siteTypeOfPctMax = $siteType1name;
					}					
				}
				$pct = $pctMax;
				$pctAll = join ";", @pct;
			}
		}
		else
		{
			$pct = $NULL_VALUE;
		}
		
		#################  Replace fields with new values  #################
		
		# Replace group type with 8mer/m8/1a conserved group type (that may not be consistent across group)
		$f[10] = $groupTypeForThisSite;
		
		# Replace species list with BL
		$f[11] = $blSelected;
		
		# Replace $groupSiteTypeSpeciesList ($f[12]) with PCT
		$f[12] = "";
		if ($pct) {	$f[12] = $pct; }

		# Add field for conservation score (1 = conserved; 0 = nonconserved)
		# based on site-type-specific BL threshold (that may not be consistent across group)
		$f[13] = $consScore;
		
		# Print all the data for this site (row of input data)
		print "$f[0]";
		for (my $i = 1; $i <= $#f; $i++)
		{
			(defined $f[$i] || $f[$i] ne $NULL_VALUE) ? print "\t$f[$i]" : print "\t";
		}
		print "\n";
	}
	
	%lineNum2data = ();
	%siteType2speciesList = ();
	@sitesTypesThisGroup = ();
}

sub lookAtBranchLength
{
	# Wrapper for getBranchLength()
	# See if we have a BL calculated for this set of species in this bin.
	# If not, we'll calculate it.

	my $refseq = $_[0];
	my $speciesList = $_[1];
	my $bl;
	
	if (! $speciesList)
	{
		print STDERR "No species list for $groupLine\n";
	}

	if ($refseq2bin{$refseq} || $refseq2bin{$refseq} eq "0")
	{
		$refseqBin = $refseq2bin{$refseq};
	}
	else
	{
		print STDERR "ERROR: No bin assignment for $refseq ; assigning to bin 1\n";
		
		# Assign 0 bin for now
		$refseqBin = 1;
	}
	
	# We already got a BL for this species list and this bin
	
	if ($speciesPlusBin2BL{"$speciesList:$refseqBin"})
	{
		$bl = $speciesPlusBin2BL{"$speciesList:$refseqBin"}
	}
	else	# We have to get the BL now
	{
		@species = split (/ /, $speciesList);
		
		if ($#species == 0)	# Only one species; BL = 0
		{
			$speciesPlusBin2BL{"$speciesList:$refseqBin"} = 0;
			
			$bl = 0;
		}
		elsif ($refseqBin || $refseqBin eq "0")
		{	
			$bl = getBranchLength($refseqBin, $speciesList);

			# Save this if we come across the same species list and bin number again 
			$speciesPlusBin2BL{"$speciesList:$refseqBin"} = $bl;
		}
		else
		{
			print STDERR "Error: No bin number for $refseq\n";
			
			$bl = $NULL_VALUE;
		}
	}
	
	return $bl;
}

sub getGroupTypeAnotherWay
{
	my @sitesTypesThisGroup = @_;
	my $siteType;
	
	my %gotThisSiteType = ();

	my @sitesTypesThisGroupUnique = grep { ! $gotThisSiteType{$_} ++ } @sitesTypesThisGroup;

	@sitesTypesThisGroupUnique = sort @sitesTypesThisGroupUnique;

	my $groupType = join "+", @sitesTypesThisGroupUnique;
	
	return $groupType;
}

sub selectGroupTypeBLforThisSite
{
	my $blAll = $_[0];
	my $blSelected;
	
	my @BLs = split /;/, $blAll;
	
	if ($BLs[0] ne $NULL_VALUE && $BLs[0] >= $siteTypeToConsThreshold{$siteType3name}) 
	{  
		$blSelected = $BLs[0];
		return "$siteType3name $blSelected";
	}
	elsif ($BLs[1] ne $NULL_VALUE && $BLs[1] >= $siteTypeToConsThreshold{$siteType2name}) 
	{  
		$blSelected = $BLs[1];
		return "$siteType2name $blSelected";
	}	
	elsif ($BLs[2] ne $NULL_VALUE && $BLs[2] >= $siteTypeToConsThreshold{$siteType1name}) 
	{  
		$blSelected = $BLs[2];
		return "$siteType1name $blSelected";
	}		
	else	# Stay with 8mer 
	{
		return "$siteType3name $BLs[0]";
	}
}

sub selectConservationForThisSite
{
	# Given the site type and BL, let's see which side of the conservation threshold we're on.

	my $bl = $_[0];
	my $siteType = $_[1];
	
	if ($siteType eq "6mer")
	{
		$cons = "";
		return $cons;
	}
	
	if (! $siteTypeToConsThreshold{$siteType})
	{
		print STDERR "No thresholds for siteType *$siteType*\n";
	}
	
	if ((! $bl && $bl != 0) || ! $siteType)
	{
		print STDERR "bl = $bl // siteType = $siteType\n";
	}
	elsif ($bl >= $siteTypeToConsThreshold{$siteType})
	{
		# It's convserved
		$cons = "x";
	}
	else
	{
		# It's poorly conserved
		$cons = "";
	}
	return $cons;
}


sub getPCT
{
	# Get the PCT from site type, miRNA family, and BL

	my $siteType = $_[0];
	my $MiRNA_ID = $_[1];
	my $branchLength = $_[2];
	
	my $pct;
	
	if ($branchLength =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
	{

		if ( $thisMiRNAhasPCT{$MiRNA_ID} && $branchLength > 0 )
		{
			# Use miRNA family, site type, and BL to calculate PCT
			# $pct = calculatePCTthisBL($MiRNA_ID, $siteType, $branchLength);
			$pct = ($siteType ne "6mer")? calculatePCTthisBL($MiRNA_ID, $siteType, $branchLength) : "0.0";
		}
		elsif ( $thisMiRNAhasPCT{$MiRNA_ID} && $branchLength == 0 )
		{
			$pct = "0.0";
		}
		else
		{
			$pct = $NULL_VALUE;
		}
	}
	else
	{
		print STDERR "Skipped trying to round *$branchLength*: Is this a problem?\n";
		$pct = "EMPTY";
	}
	return $pct;
}

sub calculatePCTthisBL
{
	# Use miRNA family, site type, and BL to calculate PCT
	
	my ($miRNAfamily, $siteType, $BL) = @_;
	
	$b0 = $familyPlusType2coeff{$miRNAfamily}{$siteType}{0};
	$b1 = $familyPlusType2coeff{$miRNAfamily}{$siteType}{1};
	$b2 = $familyPlusType2coeff{$miRNAfamily}{$siteType}{2};
	$b3 = $familyPlusType2coeff{$miRNAfamily}{$siteType}{3};
	
	if (! $b0 || !$b1 || !$b2 || !$b3)
	{
		print STDERR "* Coeffs for $miRNAfamily (type=$siteType): $b0, $b1, $b2, $b3 *\n";
	}
	
	my $pct = $b0 + ( $b1 / (1 + $eConstant ** ( (0 - $b2) * $BL + $b3)));
	
	$pct = sprintf ("%.4f", $pct);
	
	# If it's negative, set to 0
	if ($pct < 0) { $pct = "0.0"; }

	return $pct;
}

sub getBranchLength
{
	# Calculate branch length from a list of species IDs
	#
	# This section of code largely written by Robin Friedman, MIT

	my $bin = $_[0];
	my $speciesList = $_[1];
	my @orgs = split / /, $speciesList;
	
	# my $treeFile = "$pathToTreesParams/tree_bin_${bin}.txt";
	if ($bin < 10) { $bin = "0${bin}"; }
	$treeFile = "$pathToTreesParams/Tree.bin_${bin}.txt";
	
	if (! $currentBin || $bin ne $currentBin)
	{
		# Read a new 3' UTR tree for this bin 
		# print STDERR "Reading bin $bin tree ($treeFile)\n";
		READ_TREE($treeFile);
	}
	
	$currentBin = $bin;
	
	my %ref_ancestors;	# keys are internal ids, values are objects
	my %ref_cumul_dist;	# keys are internal ids, values 

	# are cumulative distance from node1 to given node
	# my $place = $nodes{$refGenome};		# start at node1
	my $place = $nodes{$orgs[0]};		# start at node1; 30 Sept 2014
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
	$total_dist = sprintf ("%.4f", $total_dist);
	
	# Return the BL
	return $total_dist;
}

sub READ_TREE 
{
	# Read a precomputed phylogenetic tree

	my $treefile = shift;
	
	$input = new Bio::TreeIO(-file   => $treefile,
					-format => "newick");
	$tree = $input->next_tree;

	foreach my $org (@all_orgs) 
	{
		($nodes{$org}) = $tree->find_node(-id => $org);
	}
}

