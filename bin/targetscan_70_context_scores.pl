#! /usr/bin/env perl
use warnings;

#######################################################################
#
# Copyright(c) 2007-2018 Whitehead Institute for Biomedical Research.
#              All Rights Reserved
#
# Author: Joe Rodriguez, Robin Ge, Kim Walker, and George W. Bell
#         Bioinformatics and Research Computing
#         wibr-bioinformatics@wi.mit.edu
#
# Comment: This program calculates TargetScan context scores.
#          It produces output as displayed in TargetScan Release 7.0
#
# This code is available from http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70
#
# Version: 5.3 (October 2010)
#          - Correct bugs in 5.0: AU score (A1 for 8mer; <30nt denominator; <30nt $utrUp)
#          - Add code for three site types
#              Nevertheless, we need this code to produce miRNA/UTR alignments
#
# Version: 6.0 (July 2011)
#          - Update to calculate context+ scores (Garcia et al.)
#
# Version: 7.0 (July 2014 - July 2015)
#          - Update to calculate context++ scores (Agarwal et al.)
#
# Version: 7.01 (25 March 2016)
#          - Update to correct errors in the calculations of the ORF 8mer and offset 6mer contributions.

# Version: 7.02 (17 August 2018)
#          - Convert mature miRNA sequence to uppercase (if it's not that way in the input file)
#
#######################################################################

# Basic ideas:
#
# 1 - Get names of all input files: 
#		- predicted targets (output from targetscan_70.pl)
#		- UTRs (used by targetscan_70.pl to predicted targets)
#		- mature miRNAs with miRNA name, family name, seed region, species ID, and mature sequence
# 2 - Get mature miRNAs and UTRs into memory
# 3 - Read through predicted targets file, and for each site 
#		- get all miRNAs in that miRNA family and species [if there are none, we can't calculate context score]
#		- extract short subsequence from UTR to use for predicted consequential pairing
#		- do and score alignment 
#		- extract long subsequence from UTR and get local AU contribution
#		- get position contribution
#		- sum 4 contributions to get total context score
#		- and print out these data
# 4 - Once all context scores have been calculated, go back and calculate context score percentiles.
# 5 - Print out all data, including percentile ranks, into final file

# Sample command
# ./targetscan_70_context_scores.pl miR_for_context_scores.sample.txt UTR_Sequences_sample.txt targetscan_70_output.BL_PCT.txt ORF_Sequences_sample.lengths.txt ORF_8mer_counts_sample.txt Targets.BL_PCT.context_scores.txt

use warnings;
# Needed for percentile rank
use POSIX qw(floor ceil);

our ($PAIRING_SCORE);

# Required!!! RNAplfold (from the ViennaRNA Package 2 -- http://www.tbi.univie.ac.at/RNA/documentation.html)

########################  Constants  ########################

# Link human-readable site types to site type code
$siteTypeDescription2num{'7mer-1a'} = 1;
$siteTypeDescription2num{'7mer-m8'} = 2;
$siteTypeDescription2num{'8mer-1a'} = 3;
$siteTypeDescription2num{'6mer'} = 4;

$siteTypeNum2description{1} = "7mer-1a";
$siteTypeNum2description{2} = "7mer-m8";
$siteTypeNum2description{3} = "8mer-1a";
$siteTypeNum2description{4} = "6mer";

# Required file of TA and SPS values (from Garcia et al., Supp Data 5)
$TA_SPS_FILE = "TA_SPS_by_seed_region.txt";

# Agarwal coefficients file
$AgarwalParamFile = "Agarwal_2015_parameters.txt";

# AIRs (UTR profiles) file
$AIRsFile = "All_cell_lines.AIRs.txt";

# Minimum distance to end of CDS.  If site is closer, context scores are not calculated
$MIN_DIST_TO_CDS = 15;
$TOO_CLOSE_TO_CDS = "too_close";
$DESIRED_UTR_ALIGNMENT_LENGTH = 23;

# Number of digits after decimal point for context scores
$DIGITS_AFTER_DECIMAL = 3;

# Print a dot every this number of input lines (so we can see progress)
$DOT_EVERY_THIS_NUM_LINES = 50;

# Maximum values for each site type
$maxContextScore{1} = -0.01;	# 7mer-1a
$maxContextScore{2} = -0.02;	# 7mer-m8
$maxContextScore{3} = -0.03;	# 8mer-1a
$maxContextScore{4} = 0;	# 6mer

# Produce all contributions to context++ score? (0=no; 1=yes)
$PRINT_CS_CONTRIBUTIONS = 1;

# Should we set all AIRs to 1 [1] or leave them as they are [0]?
$setAIRs_to_1 = 0;
# AIRs of 0 (even due to rounding errors) cause problems so set a min value
$minAIR = 0.0001;

# Directory with RNAplfold input and output
# Needs to contain one dir per species ID, with RNAplfold for that species' UTRs inside
$RNAplfold_IN_OUT = "RNAplfold_in_out";

# Reference species (of which affected isoform ratios (AIRs) were determined
$REF_SPECIES = 9606;

# Only calculate context scores for these species
@SPECIES = qw(10090 10116 13616 8364 9031 9544 9598 9606 9615 9913);

########################  Beginning of real code  ########################

getUsage();
getFileFormats();
checkArguments();

# Read coefficients from Agarwal et al. Supp Table 3
readAgarwalParameters($AgarwalParamFile);

print STDERR "Using $REF_SPECIES as the reference species (as set by \$REF_SPECIES).\n";
foreach $species (@SPECIES)
{
	$useSpecies{$species} = 1;
}

# Get mature miRNA data
readMiRNAs();

# Read file of aligned UTRs
our %utrLengthScalingFactor;
readUTRs();

# Get AIRs (UTR profiles)
readIsoformRatios($AIRsFile);

# Run RNAplfold for all UTRs
runRNAplfold_all_UTRs();

# Read files of other data
my %orf2length;
readORFdataFiles($orfLengthsFile, $orf8mersFile);

# Get mature miRNA data
read_TA_SPS($TA_SPS_FILE);

# Get a set of empty strings needed for correction of pairing
getReplaceLength(20);

# Open file for output
open (CONTEXT_SCORES_OUTPUT_TEMP, ">$contextScoresOutputTemp") || die "Cannot open $contextScoresOutputTemp for writing: $!";

if ($PRINT_CS_CONTRIBUTIONS)
{
	$contextScoreFieldNum = 27;	# 0-based
	$percentileRankFieldNum = 28;	# 0-based
	$weightedContextScoreFieldNum = 30;	# 0-based
	$weightedPercentileRankFieldNum = 31;	# 0-based
}
else
{
	$contextScoreFieldNum = 6;	# 0-based
	$percentileRankFieldNum = 7;	# 0-based
	$weightedContextScoreFieldNum = 9;	# 0-based
	$weightedPercentileRankFieldNum = 10;	# 0-based
}

# Read file of targets predicted by targetscan_70.pl
# and get context score(s) for sites line by line
readTargets();

# Read all context scores we calculated and get percentile ranks
readContextScoresFromFiles("$contextScoresOutputTemp");
getPercentileRanks();
# Read all data again, add percentile ranks, and print data to final file
addPercentileRanksToOtherData($contextScoresOutputTemp, $contextScoreFileOutput);

# Delete old RNAplfold files
unlink $contextScoresOutputTemp ;
print STDERR "Deleted temporary file ($contextScoresOutputTemp)\n";

# Report again if we're missing RNAplfold files
@missingRNAplfoldFiles = keys %missingRNAplfoldFile;
$numMissingRNAplfoldFiles = $#missingRNAplfoldFiles + 1;
if ($numMissingRNAplfoldFiles)
{
	print STDERR "\n! Note that we couldn't run RNAplfold (or find the output files) for $numMissingRNAplfoldFiles sequences,\nso these will have no SA contributions!\n";
}

print STDERR "\nAll done!  -  See $contextScoreFileOutput for output.\n\n";

########################  Subroutines  ########################

sub read_TA_SPS
{
	my $TA_SPS_file = $_[0];
	open (TA_SPS_FILE, $TA_SPS_file) || die "Cannot open file of TA and SPS values by seed ($TA_SPS_file): $!";
	while (<TA_SPS_FILE>)
	{
		# Seed region	SPS (8mer and 7mer-m8)	SPS (7mer-1a)	TA
		# GAGGUAG	-9.25	-6.72	3.393

		chomp;
		if ($. > 1)	# Skip header line
		{
			my ($seedRegion, $SPS_1, $SPS_2, $TA) = split (/\t/, $_);
			
			$garcia{3}{$seedRegion}{"SPS"} = $SPS_1;
			$garcia{2}{$seedRegion}{"SPS"} = $SPS_1;
			$garcia{1}{$seedRegion}{"SPS"} = $SPS_2;
			$garcia{4}{$seedRegion}{"SPS"} = $SPS_2;
			$garcia{$seedRegion}{"TA"} = $TA;
		}
	}
}

sub readTargets
{
	# Read file with predicted targets and their BLs and PCTs, so
	#   two previous steps are required: target prediction, BL+PCT calculation
	
	my $lineNum = 0;
	
	print STDERR "Reading and processing targets file ($predictedTargetsFile)...";
	
	open (PREDICTED_TARGETS, $predictedTargetsFile) || die "Cannot open $predictedTargetsFile for reading: $!";
	while (<PREDICTED_TARGETS>)
	{
		# Gene_ID	miRNA_family_ID	species_ID	MSA_start	MSA_end	UTR_start	UTR_end	Group_num	Site_type	miRNA in this species	Group_type	Branch length	Pct	Conserved
		# ENST00000002829	CACAGUG	9913	1380	1399	674	681	18852	8mer-1a	x	8mer-1a	1.23658845	0.4668	x
		# ENST00000002829	CACAGUG	9925	1380	1399	670	677	18852	8mer-1a	 	8mer-1a	1.23658845	0.4668	x
		# ENST00000001008	GAGGUAG	13616	49	76	31	37	7812	7mer-m8	x	7mer-m8	1.2446078	0.4806	
		# ENST00000000233	CAGCAGG	10029	19	25	18	24	958	7mer-m8	 	7mer-m8	5.14948759	0.4835	x
	
		chomp;
		my @f = split (/\t/, $_);

		my $transcriptID = $f[0];
		my $miRNA_familyID = $f[1];
		my $speciesID = $f[2];
		
		if ($useSpecies{$speciesID})
		{
			my $utrStart = $f[5];	# UTR start
			my $utrEnd = $f[6];	# UTR end
			my $groupNum = $f[7];
			my $pct = $f[12]; 

			# print "$transcriptID :: $miRNA_familyID :: $f[8]\n";

			# Type of individual site (not always the same as the group type)
			my $siteType = $siteTypeDescription2num{$f[8]};

			my $familySpecies = "$miRNA_familyID\t$speciesID";

			if ( $mature_seq{$familySpecies} )
			{
				###
				###  Get AIR (affected isoform ratio for this site)
				###

				$air = getAirThisSite($transcriptID, $speciesID, $utrEnd);

				##  If there are annotated miRNAs in this miRNA family in this species,
				##  extract piece of UTR that we'll need to align with the mature miRNA

				my $subseqForAlignment = extractSubseqForAlignment($transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType);

				###
				###  Extract long subsequence from UTR and get local AU contribution
				###

				my $localAUcontribution = getLocalAU_contribution($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd);

				###
				###  Get 3' UTR length contribution (new for TargetScan 8)
				###

				$len3UTR_contribution = get_len3UTR_weighted_contribution($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd, $air);

				###
				###  Get min_dist (shortest distance to stop codon or polyA site) contribution (new for TargetScan 8)
				###

				$minDist_contribution = getMinDist_weighted_contribution($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd);

				###
				###  Get site accessibility use ViennaRNA's RNAplfold (new for TargetScan 8)
				###

				$SA_contribution = getSA_contribution($transcriptID, $speciesID, $utrStart, $siteType);

				###
				###  Get the length of this gene's ORF (new for TargetScan 8)
				###

				$orfLength_contribution = getORFlength_contribution($transcriptID, $speciesID, $siteType);

				###
				###  Get the number of 8mers in this gene's ORF (new for TargetScan 8)
				###

				$orf8mer_contribution = getORF8mer_contribution($transcriptID, $miRNA_familyID, $speciesID, $siteType);

				###
				###  Get the number of offset 6mers in this gene's UTR (new for TargetScan 8)
				###

				$offset6mer_contribution = getOffset6mer_weighted_contribution($transcriptID, $miRNA_familyID, $speciesID, $siteType, $utrEnd);

				###
				###  Get the PCT contribution (for highly conserved miRNA families) (new for TargetScan 8)
				###    so PCT needs to have been calculated in a previous step
				###

				if ($pct eq "NA")
				{ $PCT_contribution = 0; }
				else
				{
					$PCT_contribution = getPCT_contribution($siteType, $pct);
				}

				##
				##  Get list of miRNAs for this miRNA family + species
				##

				##
				##  Get all miRNAs (if any) for this family and this species			
				##

				for (my $i = 0; $i < $#{$mature_seq{$familySpecies}} + 1; $i++)
				{
					my $thisMiRNA = @{$mature_seq{$familySpecies}}[$i];

					# print "$familySpecies ==> $thisMiRNA + UTR ($subseqForAlignment)\n";

					my ($matureMiRNAid, $matureMiRNA) = split /\t/, $thisMiRNA;

					# Modify $subseqForAlignment based on length of mature miRNA
					my ($finalSubseqForAlignment, $matureMiRNAForAlignment) = 
						modifySubseqForAlignment($matureMiRNA, $matureMiRNAid, $subseqForAlignment, $siteType);

					###
					###  Get sRNA1 and sRNA8 contributions (new for TargetScan 8)
					###

					($sRNA1A_contribution, $sRNA1C_contribution, $sRNA1G_contribution, $sRNA8A_contribution, $sRNA8C_contribution, $sRNA8G_contribution) = get_sRNA1_8_contributions($matureMiRNA, $siteType);

					###
					###  Get Site8 contributions (new for TargetScan 8)
					###

					if ($siteType == 1 || $siteType == 4)	# Only relevant for 7mer-1a and 6mer sites
					{
						($site8A_contribution, $site8C_contribution, $site8G_contribution) = getSite8_contribution($subseqForAlignment, $siteType);
					}
					else
					{
						$site8A_contribution = $site8C_contribution = $site8G_contribution = 0;
					}

					###
					###  Predict consequential pairing (alignment) and get 3' pairing contribution
					###

					my ($threePrimePairing_contribution, $alignedUTR, $alignmentBars, $alignedMatureMiRNA, $threePrimePairingScore) = get3primePairingContribution($siteType, $finalSubseqForAlignment, $matureMiRNAForAlignment);

					# This site is too close to the CDS: don't take context scores for real
					if ( $utrStart < $MIN_DIST_TO_CDS )
					{
						$threePrimePairing_contribution = 
						$localAUcontribution = 
						$totalcontextScore = 
						$TA_contribution =
						$SPS_contribution = 
						$TOO_CLOSE_TO_CDS;
					}
					else
					{
						###
						###  Sum contributions to get total context score
						###

						# Get seed region for TA and SPS
						$seedRegion = substr($matureMiRNA, 1, 7);

						if ($garcia{$seedRegion}{"TA"} && $garcia{$siteType}{$seedRegion}{"SPS"})
						{
							$TA_contribution = getAgarwalContribution($siteType, "TA_3UTR", $garcia{$seedRegion}{"TA"});
							$SPS_contribution = getAgarwalContribution($siteType, "SPS", $garcia{$siteType}{$seedRegion}{"SPS"});
						}
						else	# Set these NAs to 0 for now
						{
							$TA_contribution = 0;
							$SPS_contribution = 0;
						}

						$totalcontextScore = 
							$siteType2siteTypeContribution{$siteType} + 
							$threePrimePairing_contribution + $localAUcontribution + $minDist_contribution + 
							$sRNA1A_contribution + $sRNA1C_contribution + $sRNA1G_contribution + 
							$sRNA8A_contribution + $sRNA8C_contribution + $sRNA8G_contribution +
							$site8A_contribution + $site8C_contribution + $site8G_contribution + 
							$len3UTR_contribution + $SA_contribution + $orfLength_contribution + $orf8mer_contribution + $offset6mer_contribution + 
							$TA_contribution + $SPS_contribution + $PCT_contribution;
						$totalcontextScore = sprintf("%.${DIGITS_AFTER_DECIMAL}f", $totalcontextScore);
						# Set ceiling of total score (new for TS7)
						if ($totalcontextScore > $maxContextScore{$siteType})
						{
							$totalcontextScore = "$maxContextScore{$siteType}";
						}
					}

					###
					###  Use context score and AIR to calculate weighted context score
					###

					$weightedTotalContextScore = getWeightedContextScore($totalcontextScore, $air);

					###  
					###  Print out results for this miRNA  
					###  

					print CONTEXT_SCORES_OUTPUT_TEMP "$transcriptID";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$speciesID";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$matureMiRNAid";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$siteTypeNum2description{$siteType}";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$utrStart";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$utrEnd";

					if ($PRINT_CS_CONTRIBUTIONS)
					{
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$siteType2siteTypeContribution{$siteType}";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$threePrimePairing_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$localAUcontribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$minDist_contribution";

						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA1A_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA1C_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA1G_contribution";				
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA8A_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA8C_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$sRNA8G_contribution";

						print CONTEXT_SCORES_OUTPUT_TEMP "\t$site8A_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$site8C_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$site8G_contribution";

						print CONTEXT_SCORES_OUTPUT_TEMP "\t$len3UTR_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$SA_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$orfLength_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$orf8mer_contribution";
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$offset6mer_contribution";

						if ($garcia{$seedRegion}{"TA"} && $garcia{$siteType}{$seedRegion}{"SPS"})
						{
							print CONTEXT_SCORES_OUTPUT_TEMP "\t$TA_contribution";
							print CONTEXT_SCORES_OUTPUT_TEMP "\t$SPS_contribution";
						}
						else
						{
							print CONTEXT_SCORES_OUTPUT_TEMP "\tNA";
							print CONTEXT_SCORES_OUTPUT_TEMP "\tNA";
						}
						print CONTEXT_SCORES_OUTPUT_TEMP "\t$PCT_contribution";
					}
					# $scoreList = "$siteType2siteTypeContribution{$siteType} $threePrimePairing_contribution $localAUcontribution $positionContribution $sRNA1A_contribution $sRNA1C_contribution $sRNA1G_contribution $sRNA8A_contribution $sRNA8C_contribution $sRNA8G_contribution $site8A_contribution $site8C_contribution $site8G_contribution $len3UTR_contribution $minDist_contribution $SA_contribution $orfLength_contribution $orf8mer_contribution $offset6mer_contribution $TA_contribution $SPS_contribution";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$totalcontextScore";
					print CONTEXT_SCORES_OUTPUT_TEMP "\tNA";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$air";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$weightedTotalContextScore";
					print CONTEXT_SCORES_OUTPUT_TEMP "\tNA";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$alignedUTR";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$alignmentBars";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$alignedMatureMiRNA";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$miRNA_familyID";
					print CONTEXT_SCORES_OUTPUT_TEMP "\t$groupNum";
					print CONTEXT_SCORES_OUTPUT_TEMP "\n";
				}
			}
		}
		else
		{
			# This is a site predicted by comparative genomics and there is no annotated miRNA in this family in this species
			# Ignore this site, as context score is irrelevant 
		}
	}
	# Print a dot every $DOT_EVERY_THIS_NUM_LINES input lines
	if ($lineNum % $DOT_EVERY_THIS_NUM_LINES == 0)
	{
		print STDERR ".";
	}
	
	close (PREDICTED_TARGETS);
	close (CONTEXT_SCORES_OUTPUT_TEMP);
	
	print STDERR " done\n";
}


sub readAgarwalParameters
{
	my $AgarwalParamFile = $_[0];
	
	open (COEFFS, "$AgarwalParamFile") || die "Cannot open $AgarwalParamFile for reading: $!";
	while (<COEFFS>)
	{
		chomp;
		if ($. > 1)	# Skip header line
		{
			my ($feature, $type3Coeff, $type2Coeff, $type1Coeff, $type4Coeff, $type3Min, $type2Min, $type1Min, $type4Min, $type3Max, $type2Max, $type1Max, $type4Max) = split /\t/, $_;
			
			if ($. == 2)	# Get intercept terms
			{
				$siteType2siteTypeContribution{3} = $type3Coeff;
				$siteType2siteTypeContribution{2} = $type2Coeff;
				$siteType2siteTypeContribution{1} = $type1Coeff;
				$siteType2siteTypeContribution{4} = $type4Coeff;
			}
			else	# Get parameters
			{
				$agarwal{3}{$feature}{"coeff"} = $type3Coeff;
				$agarwal{2}{$feature}{"coeff"} = $type2Coeff;
				$agarwal{1}{$feature}{"coeff"} = $type1Coeff;
				$agarwal{4}{$feature}{"coeff"} = $type4Coeff;
				
				$agarwal{3}{$feature}{"min"} = $type3Min;
				$agarwal{2}{$feature}{"min"} = $type2Min;
				$agarwal{1}{$feature}{"min"} = $type1Min;
				$agarwal{4}{$feature}{"min"} = $type4Min;

				$agarwal{3}{$feature}{"max"} = $type3Max;
				$agarwal{2}{$feature}{"max"} = $type2Max;
				$agarwal{1}{$feature}{"max"} = $type1Max;
				$agarwal{4}{$feature}{"max"} = $type4Max;
			}
		}
	}
}


sub readUTRs
{
	# Read file with info about mature miRNAs, and get this in memory
	
	print STDERR "Reading UTRs file ($UTRfile)...";
	
	open (ALIGNED_UTRS, $UTRfile) || die "Cannot open $UTRfile for reading: $!";
	while (<ALIGNED_UTRS>)
	{
		# BMP8B	9606	GUCCACCCGCCCGGC
		# BMP8B	9615	-GUG--CUGCCCACC
	
		chomp;
		my @f = split (/\t/, $_);
		
		my $transcriptID = $f[0];
		my $speciesID = $f[1];
		$haveUTRsthisSpecies{$speciesID} = 1;
		if(! $useSpecies{$speciesID}) { next; } 	# Only get species that we want
		my $seq = $f[2];
		
		# Remove gaps from alignment
		$seq =~ s/-//g;
		
		if ($seq)	# Skip any sequences that are all gaps
		{
			# Convert T's to U's
			$seq =~ s/T/U/g;
			$seq =~ s/t/u/g;

			$utr_seq{$transcriptID}{$speciesID} = $seq;
		}
		
		if ($. % 1000 == 0)
		{
			print STDERR ".";
		}
	}
	close (ALIGNED_UTRS);
	
	# Get scaling factor for each UTR (its length vs that of the ref species)
	foreach $transcriptID ( keys %utr_seq ) 
	{
		foreach $speciesID ( keys %{ $utr_seq{$transcriptID}}) 
		{
			if (length($utr_seq{$transcriptID}{$speciesID}))
			{
				$utrLengthScalingFactor{$transcriptID}{$speciesID} = length($utr_seq{$transcriptID}{$REF_SPECIES}) / length($utr_seq{$transcriptID}{$speciesID});
				# print STDERR "utrLengthScalingFactor for *$transcriptID $speciesID* => $utrLengthScalingFactor{$transcriptID}{$speciesID}\n";
			}
		}
	}
	
	print STDERR " done\n";
}

sub readMiRNAs
{
	# Read file with info about mature miRNAs, and get this in memory

	print STDERR "Reading miRNA file ($miRNAfile)...\n";
	
	open (MATURE_MIRNAS, $miRNAfile) || die "Cannot open $miRNAfile for reading: $!";
	while (<MATURE_MIRNAS>)
	{
		# let-7/98	9606	hsa-let-7a	UGAGGUAGUAGGUUGUAUAGUU
		# miR-1/206	10090	mmu-miR-1	UGGAAUGUAAAGAAGUAUGUAU
	
		chomp;
		my @f = split (/\t/, $_);
		
		# Make key of miRNA family seed region and species ID
		my $miRNA_familyID = $f[0];
		my $speciesID = $f[1];
		my $familySpecies = "$miRNA_familyID\t$speciesID";
		
		# Make hash linking $familySpecies to mature miRNA sequence + mature miRNA name 
		my $matureMiRNAname = $f[2];
		# Aug 17, 2018: Force mature miRNA to be uppercase
		my $matureMiRNAseq = uc($f[3]);
		my $matureMiRNAnameSeq = "$matureMiRNAname\t$matureMiRNAseq";
		
		# Get seed region for family
		$seedRegion = substr($matureMiRNAseq, 1, 7);
		# my $sRNA1_nt = substr($matureMiRNAseq, 0, 1);
		# my $sRNA8_nt = substr($matureMiRNAseq, 7, 1);
		# $miRNA_familyID_to_sRNA1{$miRNA_familyID} = $sRNA1_nt;
		# $miRNA_familyID_to_sRNA8{$miRNA_familyID} = $sRNA8_nt;

		# Assign family => seed region (if we haven't seen it before)
		if (! $miRNA_familyID_to_seedRegion{$miRNA_familyID})
		{
			$miRNA_familyID_to_seedRegion{$miRNA_familyID} = $seedRegion;
		}
		# Check that family has only 1 seed region
		elsif ($seedRegion ne $miRNA_familyID_to_seedRegion{$miRNA_familyID})
		{
			print STDERR "ERROR: miRNA family $miRNA_familyID seems to have more than 1 seed region: $miRNA_familyID_to_seedRegion{$miRNA_familyID}, $seedRegion\n";
			print STDERR "Please correct family definitions and re-run analysis\n";
			exit;
		}

		if ( $mature_seq{$familySpecies} )
		{
			push (@{$mature_seq{$familySpecies}}, $matureMiRNAnameSeq);
		}
		else
		{
			@{$mature_seq{$familySpecies}} = $matureMiRNAnameSeq;
		}
	}
	close (MATURE_MIRNAS);
}

sub readIsoformRatios
{
	my $AIRsFile = shift;
	open (ISOFORM_RATIOS, $AIRsFile) || die "Cannot open $AIRsFile for reading: $!";
	print STDERR "Reading affected isoform ratios (AIRs) file ($AIRsFile)...\n";
	while (<ISOFORM_RATIOS>)
	{
		chomp;
		my ($utr, $start, $end, $air) = split /\t/, $_;
		my $scaledEnd;
		if ($setAIRs_to_1) { $air = 100; }
		# Scale these ratios according to UTR lengths in different species
		foreach $speciesID (keys %haveUTRsthisSpecies)
		{
			# Look only at UTRs in our current analysis set
			if (defined $utrLengthScalingFactor{$utr}{$speciesID})
			{
				# Need to divide by scaling factor
				$scaledEnd = sprintf("%.0f", $end / $utrLengthScalingFactor{$utr}{$speciesID});
				if ($start > 1)
				{
					$scaledStart = sprintf("%.0f", $start / $utrLengthScalingFactor{$utr}{$speciesID});
				}
				else
				{
					$scaledStart = 1;
				}
				$utrScaledEndToAIR{$utr}{$speciesID}{$scaledEnd} = $air;
				$utrScaledEndToStart{$utr}{$speciesID}{$scaledEnd} = $scaledStart;
			}
		}
	}
}

sub readORFdataFiles
{
	# Get names of other 3 files and then read them
	my ($orfLengthsFile, $orf8mersFile) = @_;
	
	open (ORF_LENGTHS, $orfLengthsFile) || die "Cannot open $orfLengthsFile for reading: $!";
	print STDERR "Reading ORF lengths file ($orfLengthsFile)...\n";
	while (<ORF_LENGTHS>)
	{
		chomp;
		my ($orf, $species, $length) = split /\t/, $_;
		if($useSpecies{$species})
		{
			$orf2length{$orf}{$species} = $length;
		}
	}

	open (ORF_8MERS, $orf8mersFile) || die "Cannot open $orf8mersFile for reading: $!";
	print STDERR "Reading ORF 8mer counts file ($orf8mersFile)...\n";
	while (<ORF_8MERS>)
	{
		chomp;
		my ($orf, $species, $miRNAfamily, $count) = split /\t/, $_;
		if($useSpecies{$species})
		{
			$orfTo8merCounts{$orf}{$species}{$miRNAfamily} = $count;
		}
	}
}

sub getAirThisSite
{
	my ($transcriptID, $speciesID, $siteEnd) = @_;
	$closestAirRegion = 10000000;
	
	if (! $utrScaledEndToAIR{$transcriptID}{$speciesID})
	{
		print STDERR "No isoform ratios (AIRs) for UTR $transcriptID of species $speciesID.  Setting AIR to 1.\n";
		$airThisRegion = 1;
	}
	else
	{
		# Go through all UTR ends in UTR profile
		foreach $AIR_end ( keys %{ $utrScaledEndToAIR{$transcriptID}{$speciesID} } ) 
		{
			$AIR_start = $utrScaledEndToStart{$transcriptID}{$speciesID}{$AIR_end};

			if ( $AIR_end >= $siteEnd && $AIR_start <= $siteEnd)	# This UTR is long enough to include this site 
			{
				# print STDERR "Site in $transcriptID ($speciesID) at $siteEnd within AIR region $AIR_start - $AIR_end (AIR=$air)\n";
				$air = $utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end};
				# Divide by 100 to change 100-point scale to 0-1 fractional scale
				return (1.0*$air/100);
			}
			else
			{
				# Need to do this to take care of rounding errors in AIR scaling
				$distToNearestEnd = min(abs($AIR_end - $siteEnd), abs($AIR_start - $siteEnd));
				if ($distToNearestEnd < $closestAirRegion)
				{
					$closestAirRegion = $distToNearestEnd;
					$airClosestRegion = $utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end};
				}
			}
		}
		if (! $airClosestRegion) { $airClosestRegion = 100; }
		$airThisRegion = 1.0*$airClosestRegion/100;
		if ($airThisRegion < $minAIR) { $airThisRegion = $minAIR; }
	}
	return $airThisRegion;
}

sub getWeightedContextScore
{
	my ($contextScore, $air) = @_;
	
	$weightedTotalContextScore = log((2**($contextScore) - 1) * $air + 1) / log(2);
	$weightedTotalContextScore = sprintf("%.4f", $weightedTotalContextScore);
	return $weightedTotalContextScore;
}

sub getUsage
{
	$usage = <<EODOCS;

	Description: Calculate context scores predicted miRNA targets
		     using TargetScan methods. 

	USAGE:
		$0 miRNA_file UTR_file PredictedTargetsBL_PCT_file ORF_lengths_file ORF_8mer_counts_file ContextScoresOutput_file

	EXAMPLE:
		$0 miR_for_context_scores.sample.txt UTR_Sequences_sample.txt targetscan_70_output.BL_PCT.txt ORF_Sequences_sample.lengths.txt ORF_8mer_counts_sample.txt Targets.BL_PCT.context_scores.txt


	Required input files:
		miRNA_file       => mature miRNA data [not used by targetscan_70.pl]
		UTR_file         => aligned UTRs (same as for targetscan_70.pl)
		PredictedTargets => output from targetscan_70_BL_PCT.pl
		ORF lengths      => length of each ORF corresponding to aligned 3\' UTRs 
		ORF 8mer counts  => number of 8mer sites in ORFs of previous file
		TA_SPS_FILE      => TA and SPS parameters (same as for targetscan_70_BL_PCT.pl) 
		                    called "TA_SPS_by_seed_region.txt"
		CS++ parameters  => Parameters for context++ score model (Agarwal et al., 2015) 
		                    called "Agarwal_2015_parameters.txt"
		UTR profiles     => Affected isoform ratios (AIRs) by 3\' UTR region
		                    called "All_cell_lines.AIRs.txt"

	Output file:
		ContextScoresOutput => Lists context scores and contributions

	For a description of input file formats, type
		$0 -h

	Authors: Joe Rodriguez, Robin Ge, Kim Walker, and George W. Bell,
	         Bioinformatics and Research Computing
	Version: 7.0 
	Copyright (c) The Whitehead Institute of Biomedical Research 

EODOCS
}

sub checkArguments
{
	# Check for input and output file arguments
	if ($ARGV[0] && $ARGV[0] eq "-h")
	{
		print STDERR "$usage";
		print STDERR "$fileFormats";
		exit (0);
	}
	elsif (! $ARGV[2])
	{
		print STDERR "$usage";
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
	elsif (! -e $ARGV[2])	# PredictedTargets file not present
	{
		print STDERR "\nI can't find the file $ARGV[2]\n";
		print STDERR "which should contain the predicted targets/BL/PCT file from targetscan_70_BL_PCT.pl.\n";
		exit;
	}
	elsif (! $ARGV[3])	# Output file not given
	{
		print STDERR "\n*** You need to supply a name for the ORF length data file ***\n";
		print STDERR "$usage";
		exit;
	}
	elsif (! $ARGV[4])	# Output file not given
	{
		print STDERR "\n*** You need to supply a name for the ORF 8mer counts file ***\n";
		print STDERR "$usage";
		exit;
	}
	elsif (! $ARGV[5])	# Output file not given
	{
		print STDERR "\n*** You need to supply a name for the context++ scores output file ***\n";
		print STDERR "$usage";
		exit (0);
	}
	
	# Get the file names
	$miRNAfile = $ARGV[0];
	$UTRfile = $ARGV[1];
	$predictedTargetsFile = $ARGV[2];
	$orfLengthsFile = $ARGV[3];
	$orf8mersFile = $ARGV[4];
	$contextScoreFileOutput = $ARGV[5];
	$contextScoresOutputTemp = "$ARGV[5].tmp";
	
	if (-e $contextScoreFileOutput)
	{
		print STDERR "Should I over-write $contextScoreFileOutput [yes/no]? ";
		$answer = <STDIN>;
		if ($answer !~ /^y/i)	{ exit; }
	}
}

sub getFileFormats
{
	$fileFormats = <<EODOCS;

	** Required input files:
	
	1 - miRNA_file    => mature miRNA information
		
		contains four fields (tab-delimited):
		Family members	Seed+m8	Species ID	miRNA_ID	Mature sequence

			a. miRNA family ID/name
			b. species ID in which this miRNA has been annotated
			c. ID for this mature miRNA
			d. sequence of this mature miRNA
			
		ex:	   
		miR-1/206	9606	hsa-miR-1	UGGAAUGUAAAGAAGUAUGUAU
		let-7/98	10090	mmu-let-7a	UGAGGUAGUAGGUUGUAUAGUU
		
	2 - UTR_file      => Aligned UTRs [same format as for targetscan_70.pl]		

		contains three fields (tab-delimited):
			a. Gene/UTR ID or name
			b. Species ID for this gene/UTR (must match ID in miRNA file)
			c. Aligned UTR or gene (with gaps from alignment)
		ex:
		BMP8B	9606	GUCCACCCGCCCGGC
		BMP8B	9615	-GUG--CUGCCCACC
		
		A gene will typically be represented on multiple adjacent lines.	

	3 - PredictedTargets_file => targets (with BL and PCT) from targetscan_70_BL_PCT.pl
	
		contains 14 fields, although some fields are ignored
		
		See targetscan_70_output.BL_PCT.txt for sample

	4 - ORF lengths (for ORFs matching 3\' UTRs in UTR_file):     file with Gene/UTR ID, species ID, length
	                                                              created by targetscan_count_8mers.pl
	
	5 - ORF 8mer counts (for ORFs matching 3\' UTRs in UTR_file): file with Gene/UTR ID, species ID, miRNA family ID/name, count of sites
	
EODOCS
}

sub extractSubseqForAlignment
{
	# Extract a subsequence of a UTR to use for prediced consequential pairing

	my ($transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType) = @_;
	my ($real_start, $real_end);
	my $startEnd = "$utrStart\t$utrEnd";

	#####
	# Step 1 - Calculate coordinates needed to extract utr_seq
	#####

	# depending on the seed_type
	# you want to do different things to the start and end coordinates
	# these start and end coordinates refer to which piece of msa_sequence you want for the RNAfold process
	
	### 4 Nov 2010 -- Add 3 new site types
	if ($siteType < 5)
	{
		$real_start = $utrStart - 16;
	}
	elsif ($siteType >= 5)
	{
		$real_start = $utrStart - 19;
	}
	else
	{ print STDERR "Unrecognized site type\n"; }
	
	if ($real_start < 0)	{ $real_start = 0; }
	
	# if ($siteType == 1 || $siteType == 2 || $siteType == 6)
	if ($siteType == 1 || $siteType == 2 || $siteType == 4)
	{
		$real_end = $utrEnd + 1;
	}
	else
	{
		$real_end = $utrEnd;
	}
	if ($real_start >= $real_end)	{ $real_start = 0; }

	#####
	# Step 2 - Extract utr_seq
	#####

	# Now that we have the right coordinates
	# lets get the length of the seq to be extracted (substr needs length)
	# use $key to pull out the right seq from %utr_seq $key = symbol_tax_id

	my $length = $real_end - $real_start;
	my $seq = $utr_seq{$transcriptID}{$speciesID};
	my $subseqForAlignment = substr($seq, $real_start, $length);

	# one last thing, we need to double check the length of the sequence
	# If the sequence doesn't equal $length ($length is how long the sequence should be)
	# that means it was either at the very end or the very beginning of the utr_seq
	# so the rest of the script runs correctly, add the correct number of N's till it is $length
	
	my $subseqForAlignmentLength = length($subseqForAlignment);

	$subseqSpacer = "";
	
	if ($DESIRED_UTR_ALIGNMENT_LENGTH > $subseqForAlignmentLength)
	{
		for($i = $subseqForAlignmentLength; $i<= $DESIRED_UTR_ALIGNMENT_LENGTH; $i++)
		{
			$subseqSpacer .= "N";
		}
		$subseqForAlignment = "$subseqSpacer$subseqForAlignment";
	}

	# Site is right at the end of the UTR
	# Remove leading Ns and add trailing Ns

	if (length($seq) < ($real_start + $length))
	{
		my $lengthDiff = $real_start + $length - length($seq);

		for($i = 0; $i < $lengthDiff; $i++)
		{
			$subseqForAlignment .= "N";
			$subseqForAlignment =~ s/^NN/  /;
		}
	}
	
	# print STDERR "*$transcriptID, $miRNA_familyID, $speciesID, $utrStart, $utrEnd, $siteType* ==> $subseqForAlignment\n";
	
	return $subseqForAlignment;
}

sub modifySubseqForAlignment
{
	# Modify $subseqForAlignment based on length of mature miRNA

	my ($matureMiRNA, $matureMiRNAid, $subseqForAlignment, $siteType) = @_;
	my $subseqForAlignmentLength = length($subseqForAlignment);

	my $spacer1Length = length($matureMiRNA) - $DESIRED_UTR_ALIGNMENT_LENGTH;
	my $spacer1 = "";
	my $spacer2 = "";
	my $spacer2Length;

	for (my $i = 0; $i < $spacer1Length; $i++)
	{
		$spacer1 .= " ";
	}

	if ($spacer1Length < 0)
	{
		$spacer2Length = -$spacer1Length;
	}
	else
	{
		$spacer2Length = 0;
	}

	for (my $i = 0; $i < $spacer2Length; $i++)
	{
		$spacer2 .= " ";
	}

	if ($DESIRED_UTR_ALIGNMENT_LENGTH > $subseqForAlignmentLength)
	{
		for(my $i = $subseqForAlignmentLength; $i <= $DESIRED_UTR_ALIGNMENT_LENGTH; $i++)
		{
			$spacer2 .= " ";
		}
	}
		############

	###  Do this adjustment for 7mer-1a _and_ 8mer sites -- 13 Oct 2010 (GB)
	if ($siteType == 1 || $siteType == 3) { chop $spacer2; }

	my $finalSubseqForAlignment = $spacer1;
	$finalSubseqForAlignment .= $subseqForAlignment;

	my $matureMiRNAForAlignment = $spacer2;
	$matureMiRNAForAlignment .= reverse $matureMiRNA;

	my @forAlignment;
	push @forAlignment, $finalSubseqForAlignment, $matureMiRNAForAlignment;
	
	return @forAlignment;
}

sub getLocalAU_contribution
{
	my ($transcriptID, $speciesID, $siteType, $utrStart, $utrEnd) = @_;

	my $totalUpScore = 0;
	my $totalDownScore = 0;
	my $utrUpStart;
	my $localAUcontributionMaxRaw = 0;

	$utrSeq = $utr_seq{$transcriptID}{$speciesID};
	$utrSeqLength = length($utrSeq);

	my $utr5 = $utrStart - 1;
	my $utr3 = $utrSeqLength - $utrEnd;
	
	# Length of UTR to extract
	my $utrSubseqLength = 30;
	
	$utrUpStart = $utrStart - 31;

	# 10 Nov 10 -- Fix short subsequences near 5' end of UTR
	if ($utrUpStart < 0)
	{
		$utrUpStart = 0;
		$utrSubseqLength = $utrStart - 1;
	}

	my $utrDownStart = $utrEnd;
	
	# Get 30 nt upstream of site
	my $utrUp = substr($utrSeq, $utrUpStart, $utrSubseqLength);
	# Get 30 nt downstream of site
	my $utrDown = substr($utrSeq, $utrDownStart, 30);

	if ( length($utrUp) < 30)
	{
		# Too close to 5' (CDS) end of UTR to get all subsequence
		# print "$transcriptID, $speciesID, $siteType, $utrStart, $utrEnd ==> too close to 5' end; getting score that we can\n";
	}
	elsif ( length($utrDown) < 30)
	{
		# Too close to 3' end of UTR to get all subsequence
		# print "$transcriptID, $speciesID, $siteType, $utrStart, $utrEnd ==> too close to 5' end; getting score that we can\n";
	}

	my @upScores = ();
	my @downScores = ();
	
	# Make an array starting with the position right before the site
	# Make it all uppercase too
	@utrUp3to5 = split(//, reverse (uc ($utrUp)));

	# Get upstream score
	for (my $i = 0; $i <= $#utrUp3to5; $i++)
	{
		if ($siteType == 2 || $siteType == 3)
		{ $scoreThisPos = 1 / ($i + 1); }
		else
		{ $scoreThisPos = 1 / ($i + 2); }
		if (($utrUp3to5[$i] eq 'U') || ($utrUp3to5[$i] eq 'A'))
		{
			$totalUpScore += $scoreThisPos;
			push @upScores, $scoreThisPos;
		}
		$localAUcontributionMaxRaw += $scoreThisPos;
	}

	@utrDown5to3 = split(//, uc ($utrDown));

	# Get downstream score
	for ($i = 0; $i <= $#utrDown5to3; $i++)
	{
		if ($siteType == 3 || $siteType == 1)
		{ $scoreThisPos = 1 / ($i + 2); }
		else
		{ $scoreThisPos = 1 / ($i + 1); }

		if (($utrDown5to3[$i] eq 'U') || ($utrDown5to3[$i] eq 'A'))
		{
			$totalDownScore += $scoreThisPos;
			push @downScores, $scoreThisPos;
		}
		$localAUcontributionMaxRaw += $scoreThisPos;
	}
	
	# Get total score and calculate regression
	$totalLocalAUscoreRaw = $totalUpScore + $totalDownScore;
	# GB 13 Oct 2010 (Prevent dividing by 0)
	if ($localAUcontributionMaxRaw != 0)
	{
		$totalLocalAUscoreFraction = sprintf("%.${DIGITS_AFTER_DECIMAL}f", $totalLocalAUscoreRaw / $localAUcontributionMaxRaw);
	}
	else
	{
		$totalLocalAUscoreFraction = 0;
	}
	
	$totalLocalAUscoreRegression = getAgarwalContribution($siteType, "Local_AU", $totalLocalAUscoreFraction);		
	return "$totalLocalAUscoreRegression";
}

sub getMinDist_weighted_contribution
{
	my ($transcriptID, $speciesID, $siteType, $siteStart, $siteEnd) = @_;	
	my @nearestEnds = ();
	my @nearestEndAirs = ();

	foreach $AIR_end ( keys %{ $utrScaledEndToAIR{$transcriptID}{$speciesID} } ) 
	{
		if ( $AIR_end >= $siteEnd)
		{
			$distTo5primeEndOfUTR = $siteStart - 1;
			$distTo3primeEndOfUTR = $AIR_end - $siteEnd;

			if ($distTo5primeEndOfUTR <= $distTo3primeEndOfUTR)
			{ $distToNearestEndOfUTR = $distTo5primeEndOfUTR; }
			else
			{ $distToNearestEndOfUTR = $distTo3primeEndOfUTR; }

			push @nearestEnds, $distToNearestEndOfUTR;
			push @nearestEndAirs, $utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end};
		}
	}
	
	$AIRs_sum = 0;
	$productSum = 0;
	for (my $i=0; $i<=$#nearestEnds; $i++)
	{
		$productSum += $nearestEnds[$i] * $nearestEndAirs[$i];
		$AIRs_sum += $nearestEndAirs[$i];
	}
	if ($AIRs_sum != 0)
	{
		$weightedMean = sprintf("%.0f", $productSum / $AIRs_sum);
	}
	else
	{
		$weightedMean = 0;
	}
	$distToNearestEndOfUTR = $weightedMean;
	
	# Transform with log10 (as of TS7)
	if (isNonzeroNumber($distToNearestEndOfUTR))	# Is it a number?
	{ $log10distToNearestEndOfUTR = log($distToNearestEndOfUTR) / log(10); }
	else
	{
		# print STDERR "distToNearestEndOfUTR from row $. is *$distToNearestEndOfUTR*, which doesn't appear to be a non-zero number.  Setting log10(distToNearestEndOfUTR) to 0\n";
		$log10distToNearestEndOfUTR = 0;	
	}
	$totalPositionRegression = getAgarwalContribution($siteType, "Min_dist", $log10distToNearestEndOfUTR);	
	return "$totalPositionRegression";
}

sub get3primePairingContribution 
{
	my ($type,$utr,$mirna) = @_;

	$utr =~ s/ //g;
	$mirna =~ s/ //g;
	$utr =~ s/\n//g;
	$mirna =~ s/\n//g;
	$utr = uc $utr;
	$mirna = uc $mirna;
	$utr =~ tr/T/U/;
	$mirna =~ tr/T/U/;
	
	###  GB - 2 Nov 2010 -- What should these values be?
	my %seedinfo = (
		'utrstart'		=>		{ 1 => 8, 2 => 8, 3 => 8, 4 => 8 },# 0 based
		'mirnastart'	=>		{ 1 => 7, 2 => 8, 3 => 8, 4 => 8 },
		'offset'		=>		{ 1 => 1, 2 => 0, 3 => 1, 4 => 2 },
		'overhang'		=>		{ 1 => 1, 2 => 0, 3 => 0, 4 => 0 },
		'seedspan'		=>		{ 1 => 6, 2 => 7, 3 => 7, 4 => 6 }
		);

	$utr = reverse $utr;
	$mirna = reverse $mirna;

	my($utrNum, $mirnaNum) = (substr($utr, $seedinfo{utrstart}{$type}),substr($mirna, $seedinfo{mirnastart}{$type})); #took part off

	my $maxscore = max(length($utrNum), length($mirnaNum));

	$utrNum =~ tr/AUCGN/12345/;
	$mirnaNum =~ tr/AUCGN/12345/;

	my @UTR = split("", $utrNum);
	my @MIRNA = split("", $mirnaNum);

	my $scorehash;
	my($prevmatch,$tempscore) = (0,0);

	for (my $offset = 0; $offset < $maxscore; $offset++) 
	{
		my $score = 0;
		my $string = "";
		my $tempstring = "";
		my $bestmatch = 0;
		for (my $i = 0; $i <= ($#MIRNA - $offset) && ($i <= $#UTR); $i++)
		{	#impose mirna limit
			if((($UTR[$i] * $MIRNA[$i + $offset]) == 2) || (($UTR[$i] * $MIRNA[$i + $offset]) == 12))
			{
				if(($i + $offset - $seedinfo{overhang}{$type}>= 4) && ($i + $offset - $seedinfo{overhang}{$type} <= 7))
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += 1;
				}
				else
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += .5;
				}
				$prevmatch++;
			}
			elsif($prevmatch >= 2)
			{
				if($tempscore == $score)
				{
					$string .= $tempstring;
				}
				elsif($tempscore > $score)
				{	# WHICH ONE DO WE TAKE IF EQUAL? don't change score, and leave both...
					$bestmatch = $prevmatch;
					$string =~ s/\|/ /g;
					$string =~ s/X/ /g;
					$string .= $tempstring;
					$score = $tempscore;
				}
				else
				{
					$tempstring =~ s/\|/ /g;
					$tempstring =~ s/X/ /g;
					$string .= $tempstring;
				}
				$string .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
			else
			{
				$tempstring =~ s/\|/ /g;
				$tempstring =~ s/X/ /g;
				$string .= $tempstring;
				$string  .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
		}
		if($prevmatch >= 2)
		{
			if($tempscore == $score)
			{
				$string .= $tempstring;
			}
			elsif($tempscore > $score)
			{
				$bestmatch = $prevmatch;
				$string =~ s/\|/ /g;
				$string =~ s/X/ /g;
				$string .= $tempstring;
				$score = $tempscore;
			}
			$tempscore = 0;
			$prevmatch = 0;
		}
		$score = $score - max(0,(($offset-2)/2));
		$string =~ s/\s([\|X])\s/   /g;
		$string =~ s/^([\|X])\s/  /g;
		$string =~ s/\s([\|X])$/  /g;
		push(@{$scorehash->{$score}}, {offset => $offset,gaploc => 'top',matchstring=>$string});

		$score = 0;
		$tempscore = 0;
		$prevmatch = 0;
		$tempstring = "";
		$bestmatch = 0;
		$string = "";
		for (my $i = 0; ($i <= ($#UTR - $offset)) && ($i <= $#MIRNA); $i++)
		{
			if(($UTR[$i + $offset] * $MIRNA[$i] == 2) || (($UTR[$i + $offset] * $MIRNA[$i]) == 12))
			{	#MATCH
				if(($i - $seedinfo{overhang}{$type}>= 4) && ($i - $seedinfo{overhang}{$type}<= 7))
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += 1;
				}
				else
				{
					$tempstring .= "|";
					if($prevmatch == 0)
					{
						$tempscore = 0;
					}
					$tempscore += .5;
				}
				$prevmatch++;
			}
			elsif($prevmatch >= 2)
			{
				if($tempscore == $score)
				{
					$string .= $tempstring;
				}
				elsif($tempscore > $score)
				{	#WHICH ONE DO WE TAKE IF EQUAL? dont change score, and leave both...
					$bestmatch = $prevmatch;
					$string =~ s/\|/ /g;
					$string =~ s/X/ /g;
					$string .= $tempstring;
					$score = $tempscore;
				}
				else
				{
					$tempstring =~ s/\|/ /g;
					$tempstring =~ s/X/ /g;
					$string .= $tempstring;
				}

				$string .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
			else
			{
				$tempstring =~ s/\|/ /g;
				$tempstring =~ s/X/ /g;
				$string .= $tempstring;
				$string  .= " ";
				$tempstring = "";
				$tempscore = 0;
				$prevmatch = 0;
			}
		}
		if($prevmatch >= 2)
		{
			if($tempscore == $score)
			{
				$string .= $tempstring;
			}
			elsif($tempscore > $score)
			{
				$bestmatch = $prevmatch;
				$string =~ s/\|/ /g;
				$string =~ s/X/ /g;
				$string .= $tempstring;
				$score = $tempscore;
			}
			$tempscore = 0;
			$prevmatch = 0;
		}
		$score = $score - max(0,(($offset-2)/2));
		$string =~ s/\s([\|X])\s/   /g;
		$string =~ s/^([\|X])\s/  /g;
		$string =~ s/\s([\|X])$/  /g;
		push(@{$scorehash->{$score}}, {offset => $offset,gaploc => 'bottom',matchstring=>$string});
		$score = 0;
		$tempscore = 0;
		$prevmatch = 0;
		$tempstring = "";
		$bestmatch = 0;
		$string = "";
	}

	#### DONE ALIGNMENT, GET BEST SCORE

	foreach my $score (sort {$b <=> $a} keys %{$scorehash})
	{
		my @outputFields;
		my($i_ret,$offset_ret);
		if(($#{$scorehash->{$score}} == 1) && ($scorehash->{$score}->[0]->{offset} == 0))
		{
			$i_ret = 1;
		}
		elsif(($#{$scorehash->{$score}} > 0))
		{
			#select one with the smallest offset, followed by shortest offset.
			#search each $i, and check offset, if previous min-i is init and greater than current offset, record min $i

			for(my $i = 0; $i <= $#{$scorehash->{$score}}; $i++)
			{
				if(defined($offset_ret))
				{
					if($scorehash->{$score}->[$i]->{offset} < $offset_ret)
					{
						$i_ret = $i;
						$offset_ret = $scorehash->{$score}->[$i]->{offset};
					}
					elsif($scorehash->{$score}->[$i]->{offset} == $offset_ret)
					{
						if(($scorehash->{$score}->[$i_ret]->{gaploc} eq "bottom") && ($scorehash->{$score}->[$i]->{gaploc} eq "bottom"))
						{
							die "ERROR Two tied scores with same offset, and gaplocation";
						}
						elsif($scorehash->{$score}->[$i]->{gaploc} eq "bottom")
						{
							$i_ret = $i;
							$offset_ret = $scorehash->{$score}->[$i]->{offset};
						}
					}
				}
				else
				{
					$i_ret = $i;
					$offset_ret = $scorehash->{$score}->[$i]->{offset};
				}
			}
		}
		else
		{
			#just print the one
			$i_ret = 0;
		}
		my $i = $i_ret;

		#adding the seedmatch check

		my($utrpre, $mirnapre) = (substr($utr, 0,$seedinfo{utrstart}{$type}),substr($mirna, 0,$seedinfo{mirnastart}{$type})); #took part off
		$utrpre =~ tr/AUCGN/12345/;
		$mirnapre =~ tr/AUCGN/12345/;
		my @UTRpre = split("",$utrpre);
		my @MIRNApre = split("",$mirnapre);
		my $matchpre = (" " x ($seedinfo{overhang}{$type} + 1));
		for (my $i = 1; $i <= $seedinfo{seedspan}{$type}; $i++)
		{
			if(($UTRpre[$i + $seedinfo{overhang}{$type}] * $MIRNApre[$i] == 2) || (($UTRpre[$i + $seedinfo{overhang}{$type}] * $MIRNApre[$i]) == 12))
			{	#MATCH
				$matchpre .= "|";
			}
			else
			{
				# !!! Something wrong: print E for the error	
				# Seed region does not show alignment -- and it always should
				$matchpre .= "E";
			}
		}
		$matchpre .= (" " x $scorehash->{$score}->[$i]->{offset});
		my $string = $matchpre . $scorehash->{$score}->[$i]->{matchstring};

		my $mirnapreoff = " " x $seedinfo{overhang}{$type};

		my $utrout = substr($utr, 0,$seedinfo{utrstart}{$type}) . substr($utr, $seedinfo{utrstart}{$type});
		my $mirnaout = $mirnapreoff . substr($mirna, 0,$seedinfo{mirnastart}{$type}) . substr($mirna, $seedinfo{mirnastart}{$type});

		my $offsetstring = "-" x $scorehash->{$score}->[$i]->{offset};
		if ($scorehash->{$score}->[$i]->{gaploc} eq "top")
		{
			$utrout = substr($utr, 0,$seedinfo{utrstart}{$type}) . "$offsetstring" . substr($utr, $seedinfo{utrstart}{$type});
		}
		elsif($scorehash->{$score}->[$i]->{gaploc} eq "bottom")
		{
			$mirnaout = $mirnapreoff . substr($mirna, 0,$seedinfo{mirnastart}{$type}) . "$offsetstring" . substr($mirna, $seedinfo{mirnastart}{$type});
		}
		my $longest = max(length($utrout),length($mirnaout));
		my $postutr =  " " x ($longest - length($utrout));
		my $postmatch = " " x ($longest - length($string));
		my $postmirna = " " x ($longest - length($mirnaout));
		$utrout .= $postutr;
		$string .= $postmatch;
		$mirnaout .= $postmirna;
		$utrout = reverse($utrout);	# UTR with gaps
		$string = reverse($string); # alignmend bars
		$mirnaout = reverse($mirnaout); # miRNA with gaps

		if ($score < 3)
		{
			# Don't show consequential pairing if this score < 3
			# This is as of 17 Oct 2008

			($utrout, $string, $mirnaout) = getConsequentialPairing($utrout, $string, $mirnaout);
		}

		$threePrimePairingRegression = getAgarwalContribution($type, "3P_score", $score);	

		# Prepare output including UTR and miRNA with gaps where appropriate
		# push @outputFields, "$threePrimePairingRegression ($score)", $utrout, $string, $mirnaout, $score;
		push @outputFields, "$threePrimePairingRegression", $utrout, $string, $mirnaout, $score;

		return @outputFields;
		last;
	}
	
	return (0,0,0);
}

sub max
{
	my @list = @_;
	my $maximum = shift @list;
	foreach (@list)
	{
		if($maximum < $_)
		{
			$maximum = $_;
		}
	}
	return $maximum;
}

sub min
{
	my @list = @_;
	my $minimum = shift @list;
	foreach (@list)
	{
		if($minimum > $_)
		{
			$minimum = $_;
		}
	}
	return $minimum;
}

sub getConsequentialPairing
{
	# Don't show consequential pairing if raw 3' pairing score < 3

	my ($utrSeq, $pairing, $mirnaSeq) = @_;
	my $newSubstring;
	my $utrSeqGapless;
	my $mirnaSeqGapless;
	
	# Working with pipes is a pain
	$pairing =~ s/\|/x/g;
	
	# If there's 1 - 2 bits of complementarity in the 3' end of the miRNA-target pair,
	# remove it, replacing it with spaces
	
	if ($pairing =~ /(\S+\s+\S+)\s+(\S+)/ || $pairing =~ /(\S+)\s+(\S+)/)
	{
		my $pairingToRemove = $1;
		my $pairingToRemoveLength = length($pairingToRemove);

		if ($replaceLength{$pairingToRemoveLength})
		{
			my $newSubstring = $replaceLength{$pairingToRemoveLength};
			
			# print "replaceLength pairingToRemoveLength ==> \"$newSubstring\"\n";		
			# print "Need to replace $pairingToRemove ($pairingToRemoveLength chars) from \"$pairing\" and replace with \"$newSubstring\"\n";

			$pairing =~ s/$pairingToRemove/$newSubstring/;
			
			###
			###  Also, remove all gaps from UTR seq and miRNA seq
			###
			
			$utrSeq = removeGapAddLeadingSpace($utrSeq);
			$mirnaSeq = removeGapAddLeadingSpace($mirnaSeq);
		}
	}
	
	# change back to pipes
	$pairing =~ s/x/|/g;

	return ($utrSeq, $pairing, $mirnaSeq);
}

sub removeGapAddLeadingSpace
{
	my $seq = $_[0];
	
	my $seqLengthPre = length($seq);
	$seq =~ s/-//g;
	$seq =~ s/\.//g;
	my $seqLengthPost = length($seq);
	
	# How many gaps were removed?
	my $seqLengthDiff = $seqLengthPre - $seqLengthPost;

	# Add leading spaces to seq if gap was removed
	if ($seqLengthDiff > 0)
	{
		my $seqGapless = "";

		for ($i = 1; $i <= $seqLengthDiff; $i++)
		{
			$seqGapless .= " ";
		}
		$seqGapless .= $seq;
		$seq = $seqGapless;
	}
	return $seq;
}

sub getAgarwalContribution
{
	my ($siteType, $contributionType, $rawScore) = @_;
	
	if ($rawScore !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/)
	{
		print STDERR "In getAgarwalContribution($siteType, $contributionType, $rawScore), score (last argument) is not a number\n";
		return 0.0000000;
	}

	# Scale score in range of 0 - 1 (for selected contributions)
	if ($contributionType =~ /^TA_3UTR$|^SPS$|^Local_AU$|^3P_score$|^SA$|^Len_ORF$|^Len_3UTR$|^Min_dist$|^PCT$/)
	{
		$scaledScore = ($rawScore - $agarwal{$siteType}{$contributionType}{"min"}) / ($agarwal{$siteType}{$contributionType}{"max"} - $agarwal{$siteType}{$contributionType}{"min"});
	}
	else	# Don't scale these
	{
		$scaledScore = $rawScore;
	}
	# Calculate this contribution
	if (! defined($agarwal{$siteType}{$contributionType}{"coeff"}))
	{
		print STDERR "PROBLEM: No coeff for $contributionType (site type = $siteType)\n";
	}
	$thisContribution = $agarwal{$siteType}{$contributionType}{"coeff"} * $scaledScore;

	# print STDERR "Agarwal $contributionType (Site type = $siteType) => raw = $rawScore => scaled = $scaledScore => contribution = $thisContribution\n";

	# Round
	return sprintf("%.${DIGITS_AFTER_DECIMAL}f", $thisContribution);
}

sub readContextScoresFromFiles
{
	# Modified for public code
	my @contextScoreFiles = @_;
	foreach $contextScoreFile (@contextScoreFiles)
	{
		open (SCORES, $contextScoreFile) || die "Cannot open $contextScoreFile: $!";
		while (<SCORES>)
		{
			chomp;
			my @f = split (/\t/, $_);

			my $miRNA = $f[2];
			my $contextScore = $f[$contextScoreFieldNum];
			my $weightedContextScore = $f[$weightedContextScoreFieldNum];

			# check that $contextScore is really a number
			if ($contextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
			{
				push @{ $miRNAtoContextScores{$miRNA} }, $contextScore;
			}
			if ($weightedContextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
			{
				push @{ $miRNAtoWeightedContextScores{$miRNA} }, $weightedContextScore;
			}
		}
		close (SCORES);
	}
}

sub getPercentileRanks
{
	foreach $miRNA (sort keys %miRNAtoContextScores)
	{
		# Sort from least negative to more negative
		my @contextScoreList = sort {$b <=> $a} @{ $miRNAtoContextScores{$miRNA} };
		my @weightedContextScoreList = sort {$b <=> $a} @{ $miRNAtoWeightedContextScores{$miRNA} };

		# print "Doing $miRNA ...\n";

		# Initialize
		my $numLower = 0;

		# Go through all the context scores in the list
		for (my $i = 0; $i <= $#contextScoreList; $i++)
		{
			my $scoreListLength = $#contextScoreList + 1;

			if ($contextScoreList[$i - 1] && $contextScoreList[$i] != $contextScoreList[$i - 1])
			{
				# Get the number of sites with a lower score
				$numLower = $i;
			}
			my $percentileRank = 100 * ( $numLower ) / $scoreListLength;

			# Get the floor so we never have a percentile rank of 100
			$percentileRank = floor($percentileRank);

			$miRNAcontextScoreToPctile{$miRNA}{$contextScoreList[$i]} = $percentileRank;

			# Weighted (modified for public code)
			if ($weightedContextScoreList[$i - 1] && $weightedContextScoreList[$i] != $weightedContextScoreList[$i - 1])
			{
				# Get the number of sites with a lower score
				$numLower = $i;
			}
			my $weightedPercentileRank = 100 * ( $numLower ) / $scoreListLength;

			# Get the floor so we never have a percentile rank of 100
			$weightedPercentileRank = floor($weightedPercentileRank);

			$miRNAcontextScoreToWeightedPctile{$miRNA}{$weightedContextScoreList[$i]} = $weightedPercentileRank;		
		}
	}
}

sub addPercentileRanksToOtherData
{
	my ($contextScoresOutputTemp, $contextScoreFileOutput) = @_;
	
	open (CONTEXT_SCORES_OUTPUT_FINAL, ">$contextScoreFileOutput") || die "Cannot open $contextScoreFileOutput for writing: $!";
	print CONTEXT_SCORES_OUTPUT_FINAL "Gene ID	Species ID	Mirbase ID	Site Type	UTR start	UTR end	";
	
	if ($PRINT_CS_CONTRIBUTIONS)
	{
		print CONTEXT_SCORES_OUTPUT_FINAL "Site type contribution	3' pairing contribution	local AU contribution	Min_dist contribution	";
		print CONTEXT_SCORES_OUTPUT_FINAL "sRNA1A contribution	sRNA1C contribution	sRNA1G contribution	sRNA8A contribution	sRNA8C contribution	sRNA8G contribution	";
		print CONTEXT_SCORES_OUTPUT_FINAL "site8A contribution	site8C contribution	site8G contribution	3'UTR length contribution	SA contribution	";
		print CONTEXT_SCORES_OUTPUT_FINAL "ORF length contribution	ORF 8mer contribution	Offset 6mer contribution	TA contribution	SPS contribution	PCT contribution	";
	}
	print CONTEXT_SCORES_OUTPUT_FINAL "context++ score	context++ score percentile	AIR	weighted context++ score	weighted context++ score percentile	UTR region	UTR-miRNA pairing	mature miRNA sequence	miRNA family	Group #\n";

	# Open the file with all the data so far
	open (CONTEXT_SCORES_OUTPUT_TEMP, $contextScoresOutputTemp) || die "Cannot open $contextScoresOutputTemp FOR READING: $!";
	while (<CONTEXT_SCORES_OUTPUT_TEMP>)
	{
		chomp;
		my @f = split (/\t/, $_);

		my $miRNA = $f[2];
		my $contextScore = $f[$contextScoreFieldNum];
		my $percentileRank = $f[$percentileRankFieldNum];
		# Below was modified for public code
		my $weightedContextScore = $f[$weightedContextScoreFieldNum];
		my $weightedPercentileRank = $f[$weightedPercentileRankFieldNum];

		# Check whether this is a miRNA for which we calculated context scores
		# And make sure it's not a "too close to ORF" site
		if ($contextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
		{
			$percentileRank = $miRNAcontextScoreToPctile{$miRNA}{$contextScore};

			$f[$percentileRankFieldNum] = $percentileRank;
		}
		
		if ($weightedContextScore =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/)
		{
			$weightedPercentileRank = $miRNAcontextScoreToWeightedPctile{$miRNA}{$weightedContextScore};

			$f[$weightedPercentileRankFieldNum] = $weightedPercentileRank;
		}

		# Add the percentile rank to the rest of the data
		my $lineWithPercentileRank = join "\t", @f;

		print CONTEXT_SCORES_OUTPUT_FINAL "$lineWithPercentileRank\n";
	}
	close(CONTEXT_SCORES_OUTPUT_TEMP);
}

sub getReplaceLength
{
	my $longestString = $_[0];
	
	my $replacementString = "";

	for (my $i = 1; $i <= $longestString; $i++)
	{
		$replacementString .= " ";
	
		$replaceLength{$i} = $replacementString;
	}
}

sub get_sRNA1_8_contributions
{
	###  Given positions 1 and 8 of the miRNA (family), return values if either of these are not "U".

	my $miRNA = $_[0];
	my $siteType = $_[1];
	
	my ($sRNA1A_contribution, $sRNA1C_contribution, $sRNA1G_contribution, $sRNA8A_contribution, $sRNA8C_contribution, $sRNA8G_contribution);

	# Get positions 1 and 8 from siRNA
	my $sRNA1_nt = substr($miRNA, 0, 1);
	my $sRNA8_nt = substr($miRNA, 7, 1);

	$isPos1{"A"} = $isPos1{"C"} = $isPos1{"G"} = 0;
	$isPos1{$sRNA1_nt} = 1;

	$isPos8{"A"} = $isPos8{"C"} = $isPos8{"G"} = 0;
	$isPos8{$sRNA8_nt} = 1;

	# Initialize each contribution to 0
	$sRNA1A_contribution = $sRNA1C_contribution = $sRNA1G_contribution = $sRNA8A_contribution = $sRNA8C_contribution = $sRNA8G_contribution = 0;

	if ($sRNA1_nt ne "U")
	{
		$sRNA1A_contribution = getAgarwalContribution($siteType, "sRNA1A", $isPos1{"A"});
		$sRNA1C_contribution = getAgarwalContribution($siteType, "sRNA1C", $isPos1{"C"});
		$sRNA1G_contribution = getAgarwalContribution($siteType, "sRNA1G", $isPos1{"G"});		
	}
	if ($sRNA8_nt ne "U")
	{
		$sRNA8A_contribution = getAgarwalContribution($siteType, "sRNA8A", $isPos8{"A"});
		$sRNA8C_contribution = getAgarwalContribution($siteType, "sRNA8C", $isPos8{"C"});
		$sRNA8G_contribution = getAgarwalContribution($siteType, "sRNA8G", $isPos8{"G"});		
	}

	return ($sRNA1A_contribution, $sRNA1C_contribution, $sRNA1G_contribution, $sRNA8A_contribution, $sRNA8C_contribution, $sRNA8G_contribution);
}

sub getSite8_contribution
{
	my $subseqForAlignment = $_[0];
	my $siteType = $_[1];
	my $sitePos8 = uc(substr($subseqForAlignment, 14, 1));

	$isSitePos8{"A"} = $isSitePos8{"C"} = $isSitePos8{"G"} = 0;
	$isSitePos8{$sitePos8} = 1;
	
	$site8A_contribution = $site8C_contribution = $site8G_contribution = 0;

	if ($sitePos8 ne "U")
	{
		$site8A_contribution = getAgarwalContribution($siteType, "Site8A", $isSitePos8{"A"});
		$site8C_contribution = getAgarwalContribution($siteType, "Site8C", $isSitePos8{"C"});
		$site8G_contribution = getAgarwalContribution($siteType, "Site8G", $isSitePos8{"G"});
	}	
	return ($site8A_contribution, $site8C_contribution, $site8G_contribution);
}

sub get_len3UTR_weighted_contribution
{
	my ($transcriptID, $speciesID, $siteType, $siteStart, $siteEnd) = @_;
	my $log10utrLength;
	my @utrLengths = ();
	my @utrAirs = ();
	my $utrLengthSum = 0;
	my $AIRs_sum = 0;
			
	foreach $AIR_end ( keys %{ $utrScaledEndToAIR{$transcriptID}{$speciesID} } ) 
	{
		if ( $AIR_end >= $siteEnd)
		{
			push @utrLengths, $AIR_end;
			if ($utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end})
			{
				push @utrAirs, $utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end};
			}
			else
			{
				print STDERR "No utrScaledEndToAIR for $transcriptID $speciesID $AIR_end [line 1895]\n";
			}
		}
	}

	for (my $i=0; $i<=$#utrLengths; $i++)
	{
		if (! $utrAirs[$i]) { print STDERR "No utrAirs $i => AIRS = @utrAirs\n"; }
		$utrLengthSum += $utrLengths[$i] * $utrAirs[$i];
		$AIRs_sum += $utrAirs[$i];
	}
	if ($AIRs_sum != 0)
	{
		$utrLength = sprintf("%.0f", $utrLengthSum / $AIRs_sum);
		# print "Weighted UTR length ($speciesID): $utrLengthSum / $AIRs_sum = $utrLength\n";
	}
	else
	{
		$utrLength = 0;
	}
	
	if (isNonzeroNumber($utrLength))	# Is it a number?
	{ $log10utrLength = log($utrLength) / log(10); }
	else
	{
		# print STDERR "UTR length from row $. is *$utrLength*, which doesn't appear to be a non-zero number.  Setting log10(utrLength) to 0\n";
		$log10utrLength = 0;	
	}

	$len3UTR_contribution = getAgarwalContribution($siteType, "Len_3UTR", $log10utrLength);
	
	return $len3UTR_contribution;
}

sub getSA_contribution
{
	my ($transcriptID, $speciesID, $utrStart, $siteType) = @_;
	
	# For 6mer and 7mer-A1 sites, we should subtract 1 nt
	# so the positioning is defined the same way as 8mer and 7mer-m8 sites
	if ($siteType == 1 || $siteType == 5)
	{ $utrStart--; }
	
	# Changed for public code
	my $RNAplfold_outfile = "$RNAplfold_IN_OUT/$transcriptID.${speciesID}_lunp";
	
	# Get the scores we want
	if (-e "$RNAplfold_outfile")
	{
		$plfold = (split /\t/, `grep -A7 -P '^$utrStart\t' $RNAplfold_outfile | tail -1 | cut -f 2-`)[13];
	}
	else	# This has been changed for public code
	{
		# print STDERR "Cannot find RNAplfold output file for $transcriptID (species $speciesID).  Using 0 for plfold score.\n";
		$missingRNAplfoldFile{$RNAplfold_outfile} = 1;
		$plfold = 0;
	}

	if (! $plfold || $plfold eq "NA")
	{
		# plfold output is not a number or something else is funky; skip this contribution
		return 0;
	}
	else
	{
		if (isNonzeroNumber($plfold))	# Is it a number?
		{ $log10_plfold = log($plfold) / log(10); }
		else
		{
			print STDERR "PLfold score from row $. is *$plfold*, which doesn't appear to be a non-zero number.  Setting log10(plfold) to 0\n";
			$log10_plfold = 0;
		}
		$SA_contribution = getAgarwalContribution($siteType, "SA", $log10_plfold);
		# print "SA_Contribution $transcriptID, $speciesID, $utrStart, $siteType => $SA_Contribution\n";
		return $SA_contribution;
	}
}

sub getORFlength_contribution
{
	my ($transcriptID, $speciesID, $siteType) = @_;
	
	# Get ORF length
	$orfLength = $orf2length{$transcriptID}{$speciesID};

	# Log10 transform
	if (isNonzeroNumber($orfLength))	# Is it a number?
	{ $log10orfLength = log($orfLength) / log(10); }
	else
	{
		# print STDERR "ORF length from row $. ($transcriptID $speciesID) is *$orfLength*, which doesn't appear to be a non-zero number.  Setting log10(orfLength) to 0\n";
		$log10orfLength = 0;
	}

	$orfLength_contribution = getAgarwalContribution($siteType, "Len_ORF", $log10orfLength);
	return $orfLength_contribution;
}

sub getORF8mer_contribution
{
	my ($transcriptID, $miRNA_familyID, $speciesID, $siteType) = @_;
	
	if ($orfTo8merCounts{$transcriptID}{$speciesID}{$miRNA_familyID_to_seedRegion{$miRNA_familyID}})
	{
		$orf8mer_count = $orfTo8merCounts{$transcriptID}{$speciesID}{$miRNA_familyID_to_seedRegion{$miRNA_familyID}};
	}
	else	# If we don't have a number for it, assume it's 0 counts
	{
		$orf8mer_count = 0;
	}
	$orf8mer_contribution = getAgarwalContribution($siteType, "ORF8m", $orf8mer_count);
	return $orf8mer_contribution;
}

sub getOffset6mer_weighted_contribution
{
	my ($transcriptID, $miRNA_familyID, $speciesID, $siteType, $siteEnd) = @_;
	my @offset6mer_counts = ();
	my $offset6mer_countsSum = 0;
	my $AIRs_sum = 0;
	my $offset6mer_countsThisUTR = 0;

	# Added for public code
	if (! $utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID})
	{
		# print STDERR "Need to get offset 6mer sites for $transcriptID $speciesID $miRNA_familyID\n";
		# Have we looked at this family with this UTR before?  If not, do so now.
		getOffset6merSites($transcriptID, $speciesID, $miRNA_familyID);
	}
	
	# Are there any offset 6mers at all?
	if ($utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID})
	{	
		@offset6merSites = split ";", $utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID};
		# print "offset6merSites => $utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID} => @offset6merSites\n";
	
		@offset6mer_counts = ();
		@utrAirs = ();
	
		# Go through all UTR ends in UTR profile
		foreach $AIR_end ( keys %{ $utrScaledEndToAIR{$transcriptID}{$speciesID} } ) 
		{
			$offset6mer_countsThisUTR = 0;
			# print "AIR_end = $AIR_end\n";
			if ( $AIR_end >= $siteEnd)	# This UTR is long enough to include this site 
			{
				# How many offset-6mer sites are in a UTR of this length?
				foreach $offset6merPos (@offset6merSites)
				{
					# print "offset6merPos = $offset6merPos\n";
					if ($offset6merPos <= $AIR_end)
					{
						$offset6mer_countsThisUTR++;
					}
					# else { print "UTR of length $AIR_end won't have an o-6mer site at $offset6merPos\n"; } 			
				}
				push @offset6mer_counts, $offset6mer_countsThisUTR;			
				push @utrAirs, $utrScaledEndToAIR{$transcriptID}{$speciesID}{$AIR_end};
			}
		}
		# print "Counts:AIRs => @offset6mer_counts @utrAirs\n";
		for (my $i=0; $i<=$#offset6mer_counts; $i++)
		{
			$offset6mer_countsSum += $offset6mer_counts[$i] * $utrAirs[$i];
			$AIRs_sum += $utrAirs[$i];
		}
		if ($AIRs_sum != 0)
		{
			$offset6mer_count = sprintf("%.0f", $offset6mer_countsSum / $AIRs_sum);
			# print "Weighted o-6mer count ($transcriptID $miRNA_familyID $speciesID site at $siteEnd): $offset6mer_countsSum / $AIRs_sum = $offset6mer_count (@offset6merSites)\n";
		}
		else
		{
			$offset6mer_count = 0;
		}
	}
	else	# If we don't have a number for it, assume it's 0 counts
	{
		$offset6mer_count = 0;
	}

	$offset6mer_contribution = getAgarwalContribution($siteType, "Off6m", $offset6mer_count);
	return $offset6mer_contribution;
}

sub getPCT_contribution
{
	my ($siteType, $pct) = @_;
	$PCT_contribution = getAgarwalContribution($siteType, "PCT", $pct);
	return $PCT_contribution
}

sub isNonzeroNumber
{
	my $maybeNumber = shift;
	if (defined($maybeNumber) && $maybeNumber =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ && $maybeNumber != 0)	# Is it a non-zero number?
	# We're OK to get a log of this
	{ return 1; }
	else	# 0 or not a number
	{ return 0; }
}

sub runRNAplfold_all_UTRs
{
	$RNAplfoldInstalled = 0;
	# Check that RNAplfold is installed
	system("RNAplfold --version");
	if ($? == -1) 
	{
		print "Failed to run RNAplfold; We'll keep going but will not calculate the SA contribution.\n";
	}
	else
	{
		$RNAplfoldInstalled = 1;
	}
	
	if ($RNAplfoldInstalled)
	{
		# check for RNAplfold directory
		if (! -e $RNAplfold_IN_OUT) { print STDERR "Creating directory $RNAplfold_IN_OUT\n"; mkdir $RNAplfold_IN_OUT; }
		
		print STDERR "Running RNAplfold on UTRs (if we didn't do so before)....\n";

		use Cwd;
		my $pwd = cwd();

		# Go to the RNAplfold dir
		chdir "$RNAplfold_IN_OUT";

		# Get all UTRs
		foreach $transcriptID (sort keys %utr_seq)
		{
			foreach $species (sort keys %{ $utr_seq{$transcriptID} })
			{
				# Skip this step is we already have the RNAplfold output file
				if (! -e "$transcriptID.${species}_lunp")
				{
					$fastaFile = "$transcriptID.$species.fa";
					open (UTR_FASTA, ">>$fastaFile") || die "Cannot open $fastaFile for writing: $!";
					print UTR_FASTA ">$transcriptID.$species\n$utr_seq{$transcriptID}{$species}\n";
					my $CMD = "RNAplfold -L 40 -W 80 -u 20 < $fastaFile";
					`$CMD`;
					$psFile = "$transcriptID.${species}_dp.ps";
					unlink $psFile;		# We don't need this
				}
			}
		}

		# Go to where we had been before
		chdir "$pwd";
	}
}

sub getOffset6merSites
{
	# Add for public version
	my ($transcriptID, $speciesID, $miRNA_familyID) = @_;
	my @offset6merSites = ();
	
	# $seedRegionRevcomp = reverse_complement($miRNA_familyID);
	$seedRegionRevcomp = reverse_complement($miRNA_familyID_to_seedRegion{$miRNA_familyID});
	# Convert any Ts to Us
	$seedRegionRevcomp =~ s/T/U/g;
	$seedRegionRevcomp =~ s/t/u/g;
	$siteToFind = substr($seedRegionRevcomp, 0, 6);

	# Make sure this is case-insensitive (i)
	while ($utr_seq{$transcriptID} =~ /$siteToFind/gi)
	{
		# $start = $-[0];
		$end = $+[0];
		push @offset6merSites, $end;
	}
	if (@offset6merSites)
	{
		$utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID} = @offset6merSites;
		return 1;
	}
	else	# No sites
	{
		$utrToOffset6merSites{$transcriptID}{$speciesID}{$miRNA_familyID} = -1;
		return -1;
	}
}

sub reverse_complement 
{
	my $dna = shift;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);
	# complement the reversed DNA sequence
	$revcomp =~ tr/ACGTUacgtu/TGCAAtgcaa/;
	return $revcomp;
}

###############
