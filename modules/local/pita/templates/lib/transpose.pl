#!/usr/bin/perl

##############################################################################
##############################################################################
##
## /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/transpose.pl
##
##############################################################################
##############################################################################
##
## Written by Josh Stuart in the lab of Stuart Kim, Stanford University.
##
##  Email address: jstuart@stanford.edu
##          Phone: (650) 725-7612
##
## Postal address: Department of Developmental Biology
##                 Beckman Center Room B314
##                 279 Campus Dr.
##                 Stanford, CA 94305
##
##       Web site: http://www.smi.stanford.edu/people/stuart
##
##############################################################################
##############################################################################
##
## Written: 00/00/02
## Updated: 00/00/02
##
##############################################################################
##############################################################################

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";

use strict;
use warnings;

# Flush output to STDOUT immediately.
$| = 1;

my @flags   = (
                  [    '-q', 'scalar',     0,     1]
                , [    '-d', 'scalar',  "\t", undef]
                , ['--file', 'scalar',   '-', undef]
	        , [   '-nd', 'scalar',  "\t", undef]
              );

my %args = %{&parseArgs(\@ARGV, \@flags)};

if(exists($args{'--help'}))
{
   print STDOUT <DATA>;
   exit(0);
}

my $verbose = not($args{'-q'});
my $delim   = $args{'-d'};
my $newDelim= $args{'-nd'};
my $file    = $args{'--file'};
my $r       = 0;
my $c       = 0;
my $ncols   = -1;
my @x;
my $filep;

$verbose and print STDERR "Reading in table from '$file'...";
open($filep, $file) or die("Could not open file '$file' for reading");
while(<$filep>)
{
  chomp;
  my @row = split($delim,$_,-1);

  for($c = 0; scalar(@row) > 0; $c++)
  {
    $x[$r][$c] = shift @row;
  }
  if($c > $ncols)
  {
    $ncols = $c;
  }
  $r++;
}
close($filep);

my $nrows = $r;

$verbose and print STDERR " done ($nrows by $ncols).\n";

$verbose and print STDERR "Transposing to $ncols by $nrows...\n";
for($c = 0; $c < $ncols; $c++)
{
  for($r = 0; $r < $nrows; $r++)
  {
    my $x = defined($x[$r][$c]) ? $x[$r][$c] : '';
    print "$x";
    if($r < $nrows - 1)
      { print $newDelim; }
  }
  print "\n";
}
$verbose and print STDERR " done.\n";

exit(0);

__DATA__
syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/transpose.pl [OPTIONS] [FILE | < FILE]

Transposes a table -- flips the rows and columns so that what
were columns in the original table become the rows and what
were rows in the original table become the columns.  The original
file is assumed to have rows delimited by newlines and columns
delimited by tabs (the column delimiter can actually be set with
the -d flag; see below).

FILE: a file containing a table with each row containing a tuple
  and each field in the tuple seperated by a delimiter.  By
  default, the delimiter is assumed to be tab.

OPTIONS are:
-d DELIM: set the delimiter for the columns to DELIM (default is
          tab).

-nd NEWDELIM set the delimiter for the columns in the output to NEWDELIM (default is tab).



