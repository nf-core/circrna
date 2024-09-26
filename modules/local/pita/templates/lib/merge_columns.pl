#!/usr/bin/perl

use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";

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

my $column_1 = get_arg("1", 0, \%args);
my $column_2 = get_arg("2", 1, \%args);
my $delimiter = get_arg("d", "", \%args);
my $reverse = get_arg("r", 0, \%args);
my $sort_by_strings = get_arg("s", 0, \%args);
my $sort_by_numbers = get_arg("sn", 0, \%args);
my $unique_merge = get_arg("u", 0, \%args);

while (<$file_ref>)
{
  chop;

  my @row = split(/\t/);

  my $first = 1;
  for (my $i = 0; $i < @row; $i++)
  {
    if ($i == $column_1)
    {
      if ($first == 1) { $first = 0; } else { print "\t"; }

      if ($unique_merge == 1 and $row[$column_2] eq $row[$column_1])
      {
	print "$row[$column_1]";
      }
      elsif ($reverse == 1)
      {
	print "$row[$column_2]$delimiter$row[$column_1]";
      }
      elsif ($sort_by_strings == 1)
      {
	if ($row[$column_1] lt $row[$column_2]) { print "$row[$column_1]$delimiter$row[$column_2]"; } else { print "$row[$column_2]$delimiter$row[$column_1]"; }
      }
      elsif ($sort_by_numbers == 1)
      {
	if ($row[$column_1] < $row[$column_2]) { print "$row[$column_1]$delimiter$row[$column_2]"; } else { print "$row[$column_2]$delimiter$row[$column_1]"; }
      }
      else
      {
	print "$row[$column_1]$delimiter$row[$column_2]";
      }
    }
    elsif ($i != $column_2)
    {
      if ($first == 1) { $first = 0; } else { print "\t"; }

      print "$row[$i]";
    }
  }

  print "\n";
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/merge_columns.pl <file>

   Merges 2 columns with a specified delimiter

   -1 <num>      First column (default: 0)
   -2 <num>      Second column (default: 1)
   -d <delim>    Delimiter (default: "")
   -r            Reverse: column 2 goes before column 1
   -s:           Print the merged columns in sorted order (sort: by strings)
   -sn:          Print the merged columns in sorted order (sort: by numbers)
   -u:           If the two columns to be merged are identical, print only one item in the resulting column

