#!/usr/bin/perl

#use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";

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

my $uniquify_columns = get_arg("c", 0, \%args);
my $force_numbers_on_all = get_arg("f", 0, \%args);
my $delim = get_arg("d", "-", \%args);
my $group_columns = get_arg("g", "0",\%args);

my %counters;

my @columns = ();
if ($uniquify_columns =~ /\d/) 
{
   for (split /,/, $uniquify_columns) 
   {
      if (/(\d+)-(\d+)/)
      { 
	 push @columns, ($1..$2);
      } 
      else 
      { 
	 push @columns, $_ ;
      }
   }
   my %unique = ();
   @columns = grep { ! $unique{$_}++ } @columns;
}


my %unique = map { $_ => (); } @columns;


while(<$file_ref>)
{
  chop;

  my @values = split(/\t/);

  if ($group_columns)
  {
      my $con_value="";
      foreach my $col (@columns)
      {
	  $con_value = $con_value . "**" . $values[$col];
	  #delete $values[$col];
      }
      #print "con value: $con_value\n";
      
      my $value = $con_value;
      if ($force_numbers_on_all) 
      { 
	  $unique{$con_value} or $unique{$con_value}++ ;
      }
      if (my $suffix = $unique{$con_value}++) 
      {
	   $unique{$con_value = sprintf "%s${delim}%i", $value, $suffix}++; 
      }

      my $length= @columns;
      my $last_val= $columns[$length-1];
      my @uniq_vals= split(/\*\*/, $con_value);
      my $len_uniq= @uniq_vals;
      $values[$last_val]= $uniq_vals[$len_uniq-1];
  }
  
  else
  {
      for my $column (@columns)
      {
	  my $value = $values[$column];
	  if ($force_numbers_on_all) 
	  { 
	      $unique{$column}{$value} or $unique{$column}{$value}++ ;
	  }
	  if (my $suffix = $unique{$column}{$value}++) 
	  {
	      $unique{$column}{$values[$column] = sprintf "%s${delim}%i", $value, $suffix}++; 
	  }
      }
  }
  print join ("\t", @values) . "\n";
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/uniquify.pl <file>

   Takes in a tab delimited file and uniquifies a column: whenever 
   encountering an entry that already appeared, appends a number
   to that column "-1", "-2", ...

   -c <num>:   Columns to uniquify (zero based, default: 0), may be several columns, e.g. -c 0,2-4,7
   -g <1/0>:   Group colums inorder to uniquify them together. Default is 0.Forexamle running: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/uniquify.pl <FILE> -c 1-3 -g 1 will treat coloumns 1-3 as one and uniquify their concatenation.   -f:         Force a number even on first occurrence
               (default: first occurrence does not get a number appended)


   -d <str>:   Use <str> as delimiter instead of the default '-' character.
