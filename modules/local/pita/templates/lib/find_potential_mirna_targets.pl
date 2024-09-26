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

my $microrna_file = get_arg("f", "", \%args);
my $lengths_str = get_arg("l", "6-8", \%args);
my $gu_wobble_str = get_arg("gu", "6;0,7;1,8;1", \%args);
my $mismatches_str = get_arg("m", "6;0,7;0,8;1", \%args);
my $loops_str = get_arg("loop", "", \%args);
my $microrna_start_pairing = get_arg("mb", 2, \%args);
my $min_target_length = get_arg("min_len", 22, \%args);
my $output_random_locations = get_arg("random", 0, \%args);
my $plm = get_arg("plm", 0, \%args);


my @lengths;
if ($lengths_str =~ /[\-]/)
{
  @lengths = split(/\-/, $lengths_str);
}
else
{
  $lengths[0] = $lengths_str;
  $lengths[1] = $lengths_str;
}
my @allowed_lengths;
for (my $i = $lengths[0]; $i <= $lengths[1]; $i++)
{
  push(@allowed_lengths, $i);
}
my $num_allowed_lengths = @allowed_lengths;
print STDERR "allowed lengths = @allowed_lengths\n";

my %allowed_gu_wobbles;
my @gu_wobble = split(/\,/, $gu_wobble_str);
for (my $i = 0; $i < @gu_wobble; $i++)
{
  my @row = split(/\;/, $gu_wobble[$i]);
  $allowed_gu_wobbles{$row[0]} = $row[1];

  print STDERR "allowed_gu_wobbles{$row[0]} = $allowed_gu_wobbles{$row[0]}\n";
}

my %allowed_mismatches;
my @mismatches = split(/\,/, $mismatches_str);
for (my $i = 0; $i < @mismatches; $i++)
{
  my @row = split(/\;/, $mismatches[$i]);
  $allowed_mismatches{$row[0]} = $row[1];

  print STDERR "allowed_mismatches{$row[0]} = $allowed_mismatches{$row[0]}\n";
}

my @allowed_loops;
my @loops = split(/\,/, $loops_str);
for (my $i = 0; $i < @loops; $i++)
{
  push(@allowed_loops, $loops[$i]);
}
print STDERR "allowed_loops = @allowed_loops\n";

my @microrna_names;
my @microrna_sequences;
open(INPUT_FILE, "<$microrna_file");
while(<INPUT_FILE>)
{
  chop;

  my @row = split(/\t/);

  push(@microrna_names, $row[0]);
  $row[1] =~ s/U/T/g;
  push(@microrna_sequences, $row[1]);
}
print STDERR "MicroRNA Names: @microrna_names\n";
print STDERR "MicroRNA Sequences: @microrna_sequences\n";

my $curr_row_num = 0;

while(<$file_ref>)
{
  chop;

  my @row = split(/\t/);

  $row[1] =~ s/U/T/g;
  my $rc_sequence = &ReverseComplement($row[1]);
  my $sequence_length = length($rc_sequence);

  #print STDERR "$row[1] --> \n";
  #print STDERR "$rc_sequence\n\n";

  $curr_row_num++;
  print STDERR $curr_row_num."\n";

  if ($output_random_locations == 0)
  {
    for (my $i = $microrna_start_pairing - 1; $i <= length($rc_sequence) - $min_target_length; $i++)
    {
	#print STDERR "i = $i\n";
    	
      my $subsequence = substr($rc_sequence, $i, $allowed_lengths[$num_allowed_lengths - 1]);

      for (my $j = 0; $j < @microrna_sequences; $j++)
      {
	my @target = &IsPotentialTarget($subsequence, $microrna_sequences[$j]);
	if ($target[0] > 0)
	{
	  for (my $k = 0; $k < @target; $k += 7)
	  {
	    my $rc_start = ($i - ($microrna_start_pairing - 1));
	    my $rc_length = $target[$k] + ($microrna_start_pairing - 1);
	    my $start = $sequence_length - 1 - $rc_start;
	    my $end = $sequence_length - 1 - ($rc_start + $rc_length - 1);
	    $start++; # 1-based coordinate correction
	    $end++;   # 1-based coordinate correction
	    my $subsequence = substr($rc_sequence, $rc_start, $rc_length);
	    print "$microrna_names[$j]\t";
	    print "$row[0]\t";
	    print "$start\t";
	    print "$end\t";
	    print "$target[$k]\t";
	    print $target[$k + 1] . "\t";
	    print $target[$k + 2] . "\t";
	    print ( ($target[$k + 3] + $target[$k + 4]) > 0 ? "1" : "0" );
	    print "\t$microrna_sequences[$j]\t";
	    print "$subsequence\t";
	    print $target[$k + 6] . "\t";
	    print $target[$k + 5] . "\n";
	  }
	}
      }
    }
  }
  else
  {
    my $random_length = length($rc_sequence) - $min_target_length;
    if ($random_length >= 0)
    {
      my $sequence_position = int(rand($random_length + 1));

      my $subsequence = substr($rc_sequence, $sequence_position, $allowed_lengths[$num_allowed_lengths - 1]);
      
      for (my $j = 0; $j < @microrna_sequences; $j++)
      {
	my $start = $sequence_length - 1 - $sequence_position + 1;
	my $end = $sequence_length - 1 - ($sequence_position + $allowed_lengths[$num_allowed_lengths - 1] - 1) + 1;
	my $length = $allowed_lengths[$num_allowed_lengths - 1];

	print "$microrna_names[$j]\t$row[0]\t$start\t$end\t$length\t0\t0\t0\t$microrna_sequences[$j]\t$subsequence\t0\t0\n";
      }
    }
  }
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub IsPotentialTarget
{
  my ($subsequence, $microrna_sequence) = @_;

  my @res;

 # print STDERR "   Subseq:$subsequence  MicroRNA:$microrna_sequence\n";

  for (my $i = $allowed_lengths[$num_allowed_lengths - 1]; $i >= $allowed_lengths[0]; $i--)
  {
    my $subseq = substr($subsequence, 0, $i);
    if (length($subseq) == $i)
    {
      my @comparison_stats = &AlignSequenceOntoSequence($subseq, $microrna_sequence, $microrna_start_pairing - 1);
      my $mismatches = $comparison_stats[0] - $comparison_stats[1];
      my $gu_wobbles = $comparison_stats[1];

      my $last_pos = $comparison_stats[2];
      if (($i == $allowed_lengths[0]) or ($last_pos == 0) or ($plm == 1))
      {
  
      if (($mismatches <= $allowed_mismatches{$i}) and ($gu_wobbles <= $allowed_gu_wobbles{$i}))
      {
	push(@res, $i);
	push(@res, $mismatches);
	push(@res, $gu_wobbles);
	push(@res, 0);
	push(@res, 0);

	my @pairings = &FindPairings($subseq, $microrna_sequence, $microrna_start_pairing - 1, -1, -1);
	push(@res, $pairings[0]);
	push(@res, $pairings[1]);

	last;
      }
      }
    }
  }

  for (my $i = 0; $i < @allowed_loops; $i++)
  {
    my $seq_length = $allowed_loops[$i];
    for (my $j = 1; $j < $seq_length - 1; $j++)
    {
      my $subseq = substr($subsequence, 0, $j) . substr($subsequence, $j + 1, $seq_length - ($j + 1));
      if (length($subseq) == $seq_length - 1)
      {
	my @comparison_stats = &AlignSequenceOntoSequence($subseq, $microrna_sequence, $microrna_start_pairing - 1);
	my $mismatches = $comparison_stats[0] - $comparison_stats[1];
	my $gu_wobbles = $comparison_stats[1];

	if (($mismatches <= $allowed_mismatches{$seq_length}) and ($gu_wobbles <= $allowed_gu_wobbles{$seq_length}))
	{
	  push(@res, $seq_length);
	  push(@res, $mismatches);
	  push(@res, $gu_wobbles);
	  push(@res, 0);
	  push(@res, 1);

	  my @pairings = &FindPairings(substr($subsequence, 0, $seq_length), $microrna_sequence, $microrna_start_pairing - 1, $j, -1);
	  push(@res, $pairings[0]);
	  push(@res, $pairings[1]);

	  #print STDERR "Found loop in target at bp " . ($j + 1) . "\n";
	}
      }

      my $subseq = substr($subsequence, 0, $seq_length);
      if (length($subseq) == $seq_length)
      {
	my $microrna_subseq = substr($microrna_sequence, 0, $j + $microrna_start_pairing - 1);
	$microrna_subseq   .= substr($microrna_sequence, $j + $microrna_start_pairing);
	my @comparison_stats = &AlignSequenceOntoSequence($subseq, $microrna_subseq, $microrna_start_pairing - 1);
	my $mismatches = $comparison_stats[0] - $comparison_stats[1];
	my $gu_wobbles = $comparison_stats[1];

	if (($mismatches <= $allowed_mismatches{$seq_length}) and ($gu_wobbles <= $allowed_gu_wobbles{$seq_length}))
	{
	  push(@res, $seq_length);
	  push(@res, $mismatches);
	  push(@res, $gu_wobbles);
	  push(@res, 1);
	  push(@res, 0);

	  my @pairings = &FindPairings($subseq, $microrna_sequence, $microrna_start_pairing - 1, -1, $j);
	  push(@res, $pairings[0]);
	  push(@res, $pairings[1]);

	  #print STDERR "Found loop in mirna at bp " . ($j + 1) . "\n";
	}
      }
    }
  }

  return @res;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub FindPairings
{
  my ($target_sequence, $mirna_sequence, $mirna_skip_start, $target_skip_bp, $mirna_skip_bp) = @_;

  my @pairings;

  my $target_zeros = "";
  my $mirna_zeros = "";
  for (my $i = 0; $i < $mirna_skip_start; $i++)
  {
    if (length($target_zeros) > 0) { $target_zeros .= ","; }
    if (length($mirna_zeros) > 0) { $mirna_zeros .= ","; }

    $target_zeros .= "0";
    $mirna_zeros .= "0";
  }

  my $mirna_index = $mirna_skip_start;
  for (my $seq_index = 0; $seq_index < length($target_sequence); $seq_index++)
  {
    if ($mirna_index - $mirna_skip_start == $mirna_skip_bp)
    {
      if (length($pairings[1]) > 0 or length($mirna_zeros) > 0) { $mirna_zeros .= ","; }
      $mirna_zeros .= "0";

      $mirna_index++;
    }

    if ($seq_index == $target_skip_bp)
    {
      if (length($pairings[0]) > 0 or length($target_zeros) > 0) { $target_zeros .= ","; }
      $target_zeros .= "0";
    }
    else
    {
      my $s1 = substr($target_sequence, $seq_index, 1);
      my $s2 = substr($mirna_sequence, $mirna_index, 1);

      my $mismatch = $s1 ne $s2 ? 1 : 0;
      #my $gu_wobble = (($s1 eq "G" and $s2 eq "A") or ($s1 eq "A" and $s2 eq "G")) ? 1 : 0;
      my $gu_wobble = (($s1 eq "C" and $s2 eq "T") or ($s1 eq "A" and $s2 eq "G")) ? 1 : 0;

      if ($mismatch == 0 or $gu_wobble == 1)
      {
	$pairings[0] .= $target_zeros;
	$pairings[1] .= $mirna_zeros;
	$target_zeros = "";
	$mirna_zeros = "";

	if (length($pairings[0]) > 0 or length($target_zeros) > 0) { $pairings[0] .= ","; }
	if (length($pairings[1]) > 0 or length($mirna_zeros) > 0) { $pairings[1] .= ","; }

	$pairings[0] .= $mirna_index + 1;
	$pairings[1] .= $seq_index + $mirna_skip_start + 1;
      }
      else
      {
	if (length($pairings[0]) > 0 or length($target_zeros) > 0) { $target_zeros .= ","; }
	if (length($pairings[1]) > 0 or length($mirna_zeros) > 0) { $mirna_zeros .= ","; }

	$target_zeros .= "0";
	$mirna_zeros .= "0";
      }

      $mirna_index++;
    }
  }

  return @pairings;
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/find_potential_mirna_targets.pl <file>

   Takes in a stab file of targets and a stab file of mirnas, 
   and outputs locations of potential start sites within the input stab 
   file, where a mirna-rna interaction may occur.

   Output format:
   <mirna> <target> <start target> <end target> <seed length> <# mismatches> <#G:U wobbles> <mir loop bp> <target loop bp> <mirna seq> <target seq> <mirna pairing> <target pairing>

   -f <str>:       The file of the MicroRNAs

   -l <num1-num2>: Search for seed lengths of num1,...,num2 to the MicroRNA (default: 6-8)

   -gu <nums>:     Lengths for which G:U wobbles are allowed and number of allowed wobbles.
                   Format of nums: <length;num G:U>,<length;num G:U>,... (default: 6;0,7;1,8;1)

   -m <nums>:      Lengths for which mismatches are allowed and number of allowed mismatches
                   Format of nums: <length;num mismatches>,<length;num mismatches>,...
                   (default: 6;0,7;0,8;1)

   -loop <nums>:   Lengths for which a single loop in either the target or the microrna is allowed
                   Format of nums: <length>,<length>,... (default: none)

   -mb <num>:      MicroRNA basepair at which pairing should start (default: 2)

   -min_len <num>: Minimum length of the target when extended beyond the seed (default: 22)

   -random:        Output a single random location from each target sequence
   
   -plm:           "Prefer Long Mismatch". If given, then ambigious seeds are reported
                   as longer seeds with a mismatch. For example, an 8-mer seed with one
                   G:U wobble at the last position will be reported as an 8-mer with a
                   G:U wobble rather than the default reporting of a complete match 7-mer.
 

