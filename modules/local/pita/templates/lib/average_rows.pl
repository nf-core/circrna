#!/usr/bin/perl

use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/format_number.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libstats.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/vector_ops.pl";

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

my $source_key_str = get_arg("k", 0, \%args);
my $print_num_averaged_rows = get_arg("n", 0, \%args);
my $skip_num = get_arg("skip", 0, \%args);
my $min = get_arg("min", "", \%args);
my $max = get_arg("max", "", \%args);
my $take_max = get_arg("take_max", "", \%args);
my $take_min = get_arg("take_min", "", \%args);
my $take_true_med = get_arg("take_true_med", "", \%args);
my $take_sum = get_arg("take_sum", "", \%args);
my $take_last = get_arg("take_last", 0, \%args);
my $take_first = get_arg("take_first", -1, \%args);
my $precision = get_arg("precision", 3, \%args);
my $std = get_arg("std", "", \%args);
my $losoe = get_arg("losoe", "", \%args);
my $list = get_arg("list", "", \%args);
my $size = get_arg("size", 0, \%args);
my $delim = get_arg("delim", "\t", \%args);
my $weights_col = get_arg("w", -1, \%args);


my @source_key = split(/\,/, $source_key_str);
my %source_key_hash;
foreach my $key (@source_key)
{
  $source_key_hash{$key} = "1";
}

my %key2id;
my %occurences_by_id;
my $key_counter = 0;
my @averaged_data;
my @averaged_counts;
my @row_counts;
my $max_columns = 0;

for (my $i = 0; $i < $skip_num; $i++) { my $line = <$file_ref>; print "$line"; }

while(<$file_ref>)
{
    chomp;

    my @row = split(/\t/);

    if (@row > $max_columns) { $max_columns = @row; }

    my $key = &GetKey(\@row);
    my $id = $key2id{$key};
    if (length($id) == 0)
    {
	$id = $key_counter++;
	$key2id{$key} = $id;
	#$averaged_data[$id][0] = $key;
    }

    if ($take_first != -1)
    {
	my $occurences = $occurences_by_id{$id};
	
	if (length($occurences) == 0 or $occurences < $take_first)
	{
	    print "$_\n";

	    $occurences_by_id{$id} = $occurences + 1;
	}
    }
    else
    {
	for (my $i = 0; $i < @row; $i++)
	{
	   if (length($row[$i]) > 0 and length($source_key_hash{$i}) == 0)
	     {
		if (length($min) > 0 and $row[$i] < $min) { $row[$i] = $min; }
		if (length($max) > 0 and $row[$i] > $max) { $row[$i] = $max; }
		
		if ($take_max == 1)
		{
		    if (length($averaged_data[$id][$i]) == 0 or $averaged_data[$id][$i] < $row[$i])
		    {
			$averaged_data[$id][$i] = $row[$i];
		    }
		}
		elsif ($take_min == 1)
		{
		    if (length($averaged_data[$id][$i]) == 0 or $averaged_data[$id][$i] > $row[$i])
		    {
			$averaged_data[$id][$i] = $row[$i];
		    }
		}
		elsif ($take_last == 1)
		{
		    $averaged_data[$id][$i] = $row[$i];
		}
		elsif ($take_sum == 1)
		{
		    if (length($averaged_data[$id][$i]) == 0)
		    {		    	
		    	$averaged_data[$id][$i] = $row[$i];
		    }
		    else
		    {
		    	$averaged_data[$id][$i] += $row[$i];
		    }
		}
		elsif ($losoe == 1)
		{
		    if (length($averaged_counts[$id][$i]) == 0)
		    {
			$averaged_counts[$id][$i] = 0;
		    }
		    $averaged_data[$id][$i][$averaged_counts[$id][$i]] = $row[$i];

		}
		elsif ($std == 1)
		{
		    if (length($averaged_counts[$id][$i]) == 0)
		    {
			$averaged_counts[$id][$i] = 0;
		    }
		    $averaged_data[$id][$i][$averaged_counts[$id][$i]] = $row[$i];
		}
		elsif ($list == 1)
		{
		    if (length($averaged_counts[$id][$i]) == 0)
		    {
			$averaged_counts[$id][$i] = 0;
		    }
		    $averaged_data[$id][$i][$averaged_counts[$id][$i]] = $row[$i];

		}
		elsif ($take_true_med == 1)
		{
		    if (length($averaged_counts[$id][$i]) == 0)
		    {
			$averaged_counts[$id][$i] = 0;
		    }
		    $averaged_data[$id][$i][$averaged_counts[$id][$i]] = $row[$i];
		}
		elsif ($weights_col > -1){
		  $averaged_data[$id][$i] += $row[$i]*$row[$weights_col];
		}
		else
		{
		    $averaged_data[$id][$i] += $row[$i];
		}

		if($weights_col > -1){
		  $averaged_counts[$id][$i]+=$row[$weights_col] ;
		}
		else{
		  $averaged_counts[$id][$i]++;
		}
		$row_counts[$id][$i]++;
	    }
	    elsif (length($source_key_hash{$i}) > 0)
	    {
	      $averaged_data[$id][$i] = $row[$i];
	    }
	}
    }
}

if ($take_first == -1)
{
    for (my $i = 0; $i < $key_counter; $i++)
    {
	if ($print_num_averaged_rows == 1)
	{
	    my $max = 0;
	    for (my $j = 1; $j < $max_columns; $j++)
	    {
		if (length($row_counts[$i][$j]) > 0 and $row_counts[$i][$j] > $max)
		{
		    $max = $row_counts[$i][$j];
		}
	    }
	    print "$max\t";
	}
	
	#print "$averaged_data[$i][0]\t";
	
	for (my $j = 0; $j < $max_columns; $j++)
	{
	  if ($j > 0) { print "\t"; }

	  if (length($source_key_hash{$j}) == 0)
	  {
	    if ($take_max == 1)
	    {
	      print &format_number($averaged_data[$i][$j], $precision);
	    }
	    elsif ($take_min == 1)
	    {
	      print &format_number($averaged_data[$i][$j], $precision);
	    }
	    elsif ($take_sum == 1)
	    {
	      print &format_number($averaged_data[$i][$j], $precision);
	    }
	    elsif ($take_last == 1)
	    {
	      print $averaged_data[$i][$j];
	    }
	    elsif ($losoe == 1)
	    {
		#for(my $k=0; $k<$averaged_counts[$i][$j]; $k++)
		#{
		#print $averaged_data[$i][$j][$k].", ";
		#}
		print LogOfSumOfExps($averaged_data[$i][$j]);
	    }
	    elsif ($std == 1)
	    {
		print &format_number(&vec_std($averaged_data[$i][$j]), $precision);
	    }
	    elsif ($list == 1)
	    {
		my @arr = @{$averaged_data[$i][$j]};
		if ( ($size == 0) || ($size > length (@arr)) )
		{
			print join($delim, @arr);
		}
		else
		{
			print join($delim, @arr[0..($size-1)]);
		}
	    }
	    elsif ($take_true_med == 1)
	    {
		print &format_number(&vec_true_median($averaged_data[$i][$j]), $precision);
	    }
	    elsif (length($averaged_counts[$i][$j]) > 0 and $averaged_counts[$i][$j]>0)
	    {
	      print &format_number($averaged_data[$i][$j] / $averaged_counts[$i][$j], $precision);
	    }
	    
	  }
	  else
	  {
	    print $averaged_data[$i][$j];
	  }
	}

	print "\n";
    }
}

sub GetKey (\@)
{
    my ($row_str) = @_;

    my @row = @{$row_str};

    my $res = $row[$source_key[0]];

    for (my $i = 1; $i < @source_key; $i++)
    {
      $res .= "\t$row[$source_key[$i]]";
    }

    return $res;
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/average_rows.pl <source file>

   Average rows in <source file> that have the same key

   -k <num>:          Column of the key (default is 0)
                      NOTE: an index of multiple keys may be specified with commas (e.g., -k 1,4,5)
   -n:                Add the number of rows that went into the averaging
   -skip <num>:       Skip num rows in the source file and just print them (default: 0)
   -min <num>:        Values below num are converted to num
   -max <num>:        Values above num are converted to num
   -precision <num>:  Precision of the numbers printed (default: 3 digits)
   -take_max:         Instead of taking average for same key, take max
   -take_min:         Instead of taking average for same key, take min
   -take_sum:         Instead of taking average for same key, take sum of values
   -take_true_med:    Instead of taking average for same key, take the true median (if number of items is even, average two middle ones)
   -take_last:        Instead of taking average for same key, take last occurence
   -take_first <num>: Instead of taking average for same key, take first <num> occurences
   -losoe:            calculate Log Of Sum Of Exps for every key
   -std:              calculate standard deviation for every key
   -list:             outputs for each key the list of values (comma seperated)
   -size <num>:       if -list is used, print only the first <num> on the list (default: print all)
   -delim <chr>:      if -list is used, print the list using <chr> as delimiter (default is tab)
   -w <num>:          calculate weighted average using column <num> as weights
