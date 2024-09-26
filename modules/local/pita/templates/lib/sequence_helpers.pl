#!/usr/bin/perl

use strict;
use POSIX;

#------------------------------------------------------------------------------------------
# SIGMA Calculation based on: http://proligo2.proligo.com/Calculation/calculation.html
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# GLOBALS
#------------------------------------------------------------------------------------------
my @ALIGNMENT_HTML_COLORS = ( "cyan", "lime", "yellow", "red", "darkred" );

my %DNA;
$DNA{"A"} = 0;
$DNA{"C"} = 1;
$DNA{"G"} = 2;
$DNA{"T"} = 3;

my %H;
$H{"AA"} = -7.9;
$H{"AC"} = -8.4;
$H{"AG"} = -7.8;
$H{"AT"} = -7.2;
$H{"CA"} = -8.5;
$H{"CC"} = -8;
$H{"CG"} = -10.6;
$H{"CT"} = -7.8;
$H{"GA"} = -8.2;
$H{"GC"} = -9.8;
$H{"GG"} = -8;
$H{"GT"} = -8.4;
$H{"TA"} = -7.2;
$H{"TC"} = -8.2;
$H{"TG"} = -8.5;
$H{"TT"} = -7.9;

my %S;
$S{"AA"} = -22.2;
$S{"AC"} = -22.4;
$S{"AG"} = -21;
$S{"AT"} = -20.4;
$S{"CA"} = -22.7;
$S{"CC"} = -19.9;
$S{"CG"} = -27.2;
$S{"CT"} = -21;
$S{"GA"} = -22.2;
$S{"GC"} = -24.4;
$S{"GG"} = -19.9;
$S{"GT"} = -22.4;
$S{"TA"} = -21.3;
$S{"TC"} = -22.2;
$S{"TG"} = -22.7;
$S{"TT"} = -22.2;

my %H_End;
$H_End{"A"} = 2.3;
$H_End{"C"} = 0.1;
$H_End{"G"} = 0.1;
$H_End{"T"} = 2.3;

my %S_End;
$S_End{"A"} = 4.1;
$S_End{"C"} = -2.8;
$S_End{"G"} = -2.8;
$S_End{"T"} = 4.1;

my $S_Symmetry = -1.4;

my %epsilon1 = ('A' => 15.4, 'C' => 7.4, 'G' => 11.5, 'T' => 8.7);
my %epsilon2 = (
		'A' => {'A' => 13.7, 'C' => 10.6, 'G' => 12.5, 'T' => 11.4},
		'C' => {'A' => 10.6, 'C' => 7.3,  'G' => 9,    'T' => 7.6},
		'G' => {'A' => 12.6, 'C' => 8.8,  'G' => 10.8, 'T' => 10},
		'T' => {'A' => 11.7, 'C' => 8.1,  'G' => 9.5,  'T' => 9.4}
               );

my %deltaS = (
	      'A' => {'A' => 24,   'C' => 17.3, 'G' => 20.8, 'T' => 23.9},
	      'C' => {'A' => 12.9, 'C' => 26.6, 'G' => 27.8, 'T' => 20.8},
	      'G' => {'A' => 13.5, 'C' => 26.7, 'G' => 26.6, 'T' => 17.3},
	      'T' => {'A' => 16.9, 'C' => 13.5, 'G' => 12.9, 'T' => 24});
my %deltaH = (
	      'A' => {'A' => 9.1, 'C' => 6.5,  'G' => 7.8,  'T' => 8.6},
	      'C' => {'A' => 5.8, 'C' => 11,   'G' => 11.9, 'T' => 7.8},
	      'G' => {'A' => 5.6, 'C' => 11.1, 'G' => 11,   'T' => 6.5},
	      'T' => {'A' => 6,   'C' => 5.6,  'G' => 5.8,  'T' => 9.1});

# For cycle restriction
  my %nimblegen_cycle_table = (
     "A" => 1,
     "C" => 2,
     "G" => 3,
     "T" => 4,
     "AA" => 4,
     "AC" => 1,
     "AG" => 2,
     "AT" => 3,
     "CA" => 3,
     "CC" => 4,
     "CG" => 1,
     "CT" => 2,
     "GA" => 2,
     "GC" => 3,
     "GG" => 4,
     "GT" => 1,
     "TA" => 1,
     "TC" => 2,
     "TG" => 3,
     "TT" => 4);

#------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------
sub Reverse
{
    my ($sequence) = @_;

    my $res = "";

    for (my $i = length($sequence) - 1; $i >= 0; $i--)
    {
	$res .= substr($sequence, $i, 1);
    }

    return $res;
}

#------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------
sub ReverseComplement
{
    my ($sequence) = @_;

    my $res = reverse($sequence);

    $res =~ tr/ACGT/TGCA/;

    return $res;
}

#------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------
sub OLD_ReverseComplement
{
    my ($sequence) = @_;

    my $res = "";

    my %complement;
    $complement{"A"} = "T";
    $complement{"C"} = "G";
    $complement{"G"} = "C";
    $complement{"T"} = "A";

    for (my $i = length($sequence) - 1; $i >= 0; $i--)
    {
	my $char = substr($sequence, $i, 1);

	my $char_complement = $complement{$char};

	if (length($char_complement) > 0) { $res .= $char_complement; }
	else { $res .= $char; }
    }

    return $res;
}

#------------------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------------------
sub AlignmentFractionToHTMLColor
{
    my ($fraction, $str) = @_;

    my $res;

    if ($fraction > 0.2)
    {
	my $color;
	if ($fraction < 0.3) { $color = $ALIGNMENT_HTML_COLORS[0]; }
	elsif ($fraction > 0.8) { $color = $ALIGNMENT_HTML_COLORS[@ALIGNMENT_HTML_COLORS]; }
	else { $color = $ALIGNMENT_HTML_COLORS[int(($fraction - 0.3) / 0.1)]; }

	#print STDERR "frac=$fraction ";
	#print STDERR ($fraction - 0.3) / 0.1;
	#print STDERR " ";
	#print STDERR (int(($fraction - 0.3) / 0.1));
	#print STDERR "\n";

	$res = "<table cellspacing=\"0\" cellpadding=\"0\"><tr><td bgcolor=\"$color\">$str</td></tr></table>";
    }
    else
    {
	$res = $str;
    }

    return $res;
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeBackgroundToMatricesRatio
{
    my ($matrix_length) = @_;

    my $min = 0.00001;
    my $max = 0.99999;
    my $res = ($min + $max) / 2;
    while ($max - $min > 0.00001)
    {
	my $value = log(1 - $res) - $matrix_length * log($res);

	if ($value > 0)
	{
	    $min = $res;
	}
	else
	{
	    $max = $res;
	}

	$res = ($min + $max) / 2;
    }

    return $res / (1.0 - $res);
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeNucleotideSetFraction ()
{
  my ($sequence, $nucleotide_set, $raw_counts) = @_;

  my $sequence_length = length($sequence);

  my $eval_str = "(\$sequence =~ tr/[$nucleotide_set]//)";
  my $fraction = eval($eval_str);

  return $raw_counts == 1 ? $fraction : ($fraction / $sequence_length);
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeGCFraction ()
{
  my ($sequence, $raw_counts) = @_;

  my $sequence_length = length($sequence);
  my $gc_fraction = ($sequence =~ tr/[GC]//);

  return $raw_counts == 1 ? $gc_fraction : ($gc_fraction / $sequence_length);
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeMeltingTemperature
{
  my ($sequence, $TMMethod, $salt_concentration, $oligo_concentration) = @_;

  my $res = "";

  my $sequence_length = length($sequence);
  if ($sequence_length > 5)
  {
    my $gc_fraction = &ComputeGCFraction($sequence);

    if ($TMMethod eq "IDTDNA")
    {
      my $local_H;
      $local_H += $H_End{substr($sequence, 0, 1)};
      $local_H += $H_End{substr($sequence, $sequence_length - 1, 1)};

      my $local_S;
      $local_S += $S_End{substr($sequence, 0, 1)};
      $local_S += $S_End{substr($sequence, $sequence_length - 1, 1)};

      for (my $i = 0; $i < $sequence_length - 1; $i++)
      {
	my $dinucleotide = substr($sequence, $i, 2);

	$local_H += $H{$dinucleotide};
	$local_S += $S{$dinucleotide};
      }

      $local_H *= 1000;

      if (($sequence_length % 2 == 0 and substr($sequence, 0, $sequence_length / 2) eq substr($sequence, $sequence_length / 2, $sequence_length / 2)) or 
	  (substr($sequence, 0, int($sequence_length / 2)) eq substr($sequence, int($sequence_length / 2) + 1, int($sequence_length / 2))))
      {
	$local_S += $S_Symmetry;
      }

      my $TM_1M_Na = $local_H / ($local_S + 1.987 * log($oligo_concentration / 2)); # Division by 2 is a hack to fit

      my $TM_final = 1 / ((1 / $TM_1M_Na) + (4.29 * $gc_fraction - 3.95) * 1e-5 * log($salt_concentration) + 9.4 * 1e-6 * log($salt_concentration) * log($salt_concentration));

      $TM_final = &format_number($TM_final - 273.15, 3);

      $res = "$TM_final";
    }
    elsif ($TMMethod eq "SIGMA")
    {
      my $seq = $sequence;
      my ($aCount, $cCount, $gCount, $tCount) = GetCounts($seq);
      my $GC = POSIX::floor(0.5 + (100 * ( $gCount + $cCount / 2) / length($seq)));
      my $MW = (252 * $aCount + 228 * $cCount + 268 * $gCount + 243 * $tCount + 61 * (length($seq) - 1));
      my $epsilon = 0;
      my @s = split(//, $seq);
      for (my $i = 0; $i < length($seq) - 1; $i++)
      {
         $epsilon += 2 * $epsilon2{$s[$i]}{$s[$i+1]} - $epsilon1{$s[$i+1]};
      }
      $epsilon += $epsilon1{$s[$#s]};
      my $quantity = 1000 / $epsilon;
      my $ug = $quantity * $MW /1000;

      my $TM;
      if (length($seq) <= 20)
      {
        $TM = (2 * ($aCount + $tCount) + 4 * ($gCount + $cCount));
      }
      else
      {
        my $sigma_dH = 0;
        my $sigma_dS = 0;
        for (my $i = 0; $i < (length($seq) - 1); $i++)
        {
          $sigma_dH += $deltaH{$s[$i]}{$s[$i+1]};
          $sigma_dS += $deltaS{$s[$i]}{$s[$i+1]};
        }
        $TM = (-1000 * $sigma_dH / (-10.8 - $sigma_dS + 1.987 * -23.5) + -273.15 + 16.6 * log($salt_concentration) / log(10));
      }

      $res = &format_number($TM, 3);
      #print "SALT=$salt_concentration\n";
      #print "GC = $gc_fraction\n";
      #print "SeqLength = $sequence_length\n";
    }
    elsif ($TMMethod eq "Nimblegen")
    {
      $res = &format_number(81.5 + 16.6 * log($salt_concentration) / log(10) + 41 * $gc_fraction - 600 / $sequence_length, 3);
      #print "SALT=$salt_concentration\n";
      #print "GC = $gc_fraction\n";
      #print "SeqLength = $sequence_length\n";
    }
    else
    {
      die ("TM Method $TMMethod not recognized\n");
    }
  }
  else
  {
    $res = "Too short";
  }

  return $res;
}

sub GetCounts
{
  my $seq = shift;
  my ($aCount,$cCount, $gCount, $tCount) = (0,0,0,0);
  my @s = split(//, $seq, -1);
  for (my $i = 0; $i < $#s; $i++)
  {
    if ($s[$i] eq "A") { $aCount++; }
    elsif ($s[$i] eq "C") { $cCount++; }
    elsif ($s[$i] eq "G") { $gCount++; }
    elsif ($s[$i] eq "T") { $tCount++; }
    else{ die ("probe contains illegal characters: $s[$i]\n"); }
  }
  return ($aCount,$cCount, $gCount, $tCount);
}

#------------------------------------------------------------------------
# @list ParseRanges (@ranges)
#------------------------------------------------------------------------
sub ParseRanges
{
   my @segments = split(",", "@_");
   my @fields                       = ();
   my ($i,$beg,$end,$inc)           = (0,-1,-1,1);
   my $seg;

   foreach $seg (@segments)
   {
      $beg = undef;
      $end = undef;
      $inc = undef;

      if($seg =~ /^(\d+)[-](\d+)$/)
      {
	$beg = $1;
	$end = $2;
        $inc = $beg <= $end ? +1 : -1;
      }
      elsif($seg =~ /^(\d+)$/)
      {
         $beg = $1;
         $end = $1;
         $inc = 1;
      }

      if(defined($beg) and defined($end) and defined($inc))
      {
	if ($inc > 0)
	{
	  for($i = $beg; $i <= $end; $i += $inc)
	  {
	    push(@fields, $i);
	  }
	}
	else
	{
	  for($i = $beg; $i >= $end; $i += $inc)
	  {
	    push(@fields, $i);
	  }
	}
      }
      else
      {
	die("\nParseRanges() error: flag expected <int><\,><int><\-><int><\,>...<\,><int>, found: @_. Exit process.\n");
      }
   }
   @fields;
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeSynthesisCycles
{
  my ($seq) = @_;

  my $res = 1;
  my $prev_id = -1;
  for (my $i = length($seq) - 1; $i >= 0; $i--)
  {
    my $id = $DNA{substr($seq, $i, 1)};

    if ($prev_id == -1)
    {
      $res += $id;
    }
    else
    {
      $res += $id > $prev_id ? ($id - $prev_id) : ($id - $prev_id + 4)
    }

    $prev_id = $id;
  }

  return $res;
}

#--------------------------------------------------------------------------------------------------------
# $DNA_sequence_verified_for_not_exceeding_probe_cycle_limit ProbeNotExceedsCycleLimit ($sequence)
#--------------------------------------------------------------------------------------------------------
sub ComputeNimblegenCycles
{
  my $seq = $_[0];
  my $current_cycle = 0;

  my $prev_letter = "";
  for (my $i = 0; $i < length($seq); $i++)
  {
     my $c = substr($seq,$i,1);
     $current_cycle += $nimblegen_cycle_table{$prev_letter.$c};
     $prev_letter = $c;
  }

  return $current_cycle;
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub AlignSequenceOntoSequence
{
  my ($seq1, $seq2, $seq2_offset) = @_;

  my @res;
  $res[0] = 0;
  $res[1] = 0;
  $res[2] = 0;

  for (my $i = 0; $i < length($seq1); $i++)
  {
    my $s1 = substr($seq1, $i, 1);
    my $s2 = substr($seq2, $i + $seq2_offset, 1);

    if ($s1 ne $s2)
    {
      $res[0]++;
      if ($i == length($seq1) - 1)
      {
      	$res[2]=1;  # Mismatch was at last position
      }
    }

    if (($s1 eq "C" and $s2 eq "T") or ($s1 eq "A" and $s2 eq "G"))
    {
      $res[1]++;
    }
  }

  #print STDERR "Aligning s1=$seq1 s2=$seq2 s2_offset=$seq2_offset Mismatches=$res[0] NonGUMismatches=" . ($res[0]-$res[1]) . " GU_Wobbles=$res[1]\n";

  return @res; # res[0]=#Mismatches, res[1]=#G:U wobbles
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub AlignSequences
{
  my ($seq1, $seq2, $max_offset, $max_mutations, $max_non_aligned) = @_;

  my $seq1_length = length($seq1);
  my $seq2_length = length($seq2);

  my @res;

  for (my $i = 0; $i <= $max_offset; $i++)
  {
    my $mutations = &ComputeMismatchesBetweenAlignedSequences($seq1, 0, $seq1_length - 1, $seq2, $i, $max_mutations, $max_non_aligned);

    if ($mutations != -1)
    {
      $res[0] = $i;
      $res[1] = $mutations;

      last;
    }
  }

  if (length($res[0]) == 0)
  {
    for (my $i = 1; $i <= $max_offset; $i++)
    {
      my $mutations = &ComputeMismatchesBetweenAlignedSequences($seq2, 0, $seq2_length - 1, $seq1, $i, $max_mutations, $max_non_aligned);
      
      if ($mutations != -1)
      {
	$res[0] = -$i;
	$res[1] = $mutations;

	last;
      }
    }
  }

  #print STDERR "AlignSequences between seq1=$seq1 and seq2=$seq2 max_offset=$max_offset max_mutations=$max_mutations max_non_aligned=$max_non_aligned res=[offset=$res[0] mutations=$res[1]]\n";

  return @res; # res[0]=offset_by_which_seq1_moves res[1]=#Mismatches  res[0]="" if not found
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
sub ComputeMismatchesBetweenAlignedSequences
{
  my ($seq1, $seq1_start, $seq1_end, $seq2, $seq2_start, $max_mutations, $max_non_aligned) = @_;

  my $res = 0;
  my $non_aligned = 0;

  for (my $i = $seq1_start; $i <= $seq1_end; $i++)
  {
    my $s1 = substr($seq1, $i, 1);
    my $s2 = substr($seq2, $i + $seq2_start, 1);

    if (length($s2) == 0)
    {
      $non_aligned++;

      if ($non_aligned > $max_non_aligned)
      {
	$res = -1;

	last;
      }
    }
    elsif ($s1 ne $s2)
    {
      $res++;

      if ($res > $max_mutations)
      {
	$res = -1;

	last;
      }
    }
  }

  #print STDERR "ComputeMismatchesBetweenAlignedSequences between seq1=$seq1($seq1_start,$seq1_end) and seq2=$seq2($seq2_start,...) max_mutations=$max_mutations max_non_aligned=$max_non_aligned res=$res\n";

  return $res;
}

#print ComputeBackgroundToMatricesRatio(157) . "\n";
#print ComputeMeltingTemperature("ACCCTTCAGCAGTTCCACACACCCTTCAGCAGTTCCACACACCCTTCAGCAGTTCCACAC", "Nimblegen", 0.2, 0) . "\n";
#print AlignSequences($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);

1
