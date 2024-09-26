#!/usr/bin/perl

use strict;
my $arg;

my $fin = \*STDIN;
my $line_to_extract = -1;
my $ignore_N = 0;
while(@ARGV)
{
  $arg = shift @ARGV;

  if ($arg eq '--help')
  {
    print STDOUT <DATA>;
    exit(0);
  }
  elsif (-f $arg)
  {
    open($fin,$arg) or die("Can't open file '$arg' for reading.");
  }
  elsif($arg eq '-l')
  {
    $line_to_extract = shift @ARGV;
  }
  elsif ($arg eq '-ignore_N')
  {
     $ignore_N = 1;
  }
  else
  {
    die("Bad argument '$arg' given.");
  }
}

my $name = "";
my $seq = "";
my $line;
while(<$fin>)
{
  chomp;
  if(/\S/)
  {
    if(/^[ ]*>/)
    {
      $line = 0;
      s/^[ ]*>[ ]*//;
      if (length($name) > 0 and length($seq) > 0)
      {
        $seq = uc ($seq);
	if ($ignore_N == 0 or index($seq, "N") == -1)
	{
	   print "$name\t$seq\n";
	}
      }
      $name = $_;
      $seq = "";
    }
    else
    {
      $line++;

      if ($line_to_extract == -1 or $line == $line_to_extract)
      {
	s/[ ]//g;
	$seq .= $_;
      }
    }
  }
}

if (length($seq) > 0)
{
  $seq = uc ($seq);
  print "$name\t$seq\n";
}

__DATA__

syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/fasta2stab.pl [OPTIONS] < FASTA

  -l <num>: Extract only line <num> of the fasta from each sequence
            (useful for parsing alignments given in fasta)

  -ignore_N:  Do not print sequences that contain N.

