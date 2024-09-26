#!/usr/bin/perl

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";


use strict;

my @cols;
my $col;
my $col_max=10000;
my @append_cols;
my @prepend_cols;


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
my $verbose = get_arg("q", 0, \%args);
if ($verbose) {$verbose=0} else {$verbose=1}
my $skip_col= get_arg("sk",0, \%args);
my $headers= get_arg("h", 0, \%args);
my $delim_in= get_arg("d", "\t", \%args);
my $delim_out=$delim_in;
if (get_arg("di", "", \%args)) {$delim_in= get_arg("di", "\t", \%args)}
if (get_arg("do", "", \%args)) {$delim_out= get_arg("do", "\t", \%args)}
my $col_ranges= get_arg("f", "1", \%args);
my $multiple_delim_ins= get_arg("m", 0, \%args);
my $suppress_blanks= get_arg("s", 0, \%args);
my $tight= get_arg("t", 0, \%args);
my $invert= get_arg("i", 0, \%args);
my $get_cols_by_header=get_arg("n","",\%args);
my $zerobased= get_arg("0", 0, \%args);
my $preserve_empty_values = get_arg("e", 0, \%args);
my $name_file = get_arg("file", "", \%args);
my $line=0;

# print STDERR join(",",@cols), "\n";
# @cols = (@prepend_cols, @cols, @append_cols);

if ($name_file ne ""){
  open (NFILE,$name_file);
  my @n;
  while(<NFILE>){
    chomp;
    push @n,$_;
  }
  close NFILE;
  $get_cols_by_header=join(",",@n);
}

my $first_line="";
if ($get_cols_by_header ne ""){
    $first_line=<$file_ref>;
    my $l=$first_line;
    chomp $l;
    my @header=split /\t/,$l,-1;
    my %heads;
    for (my $i=0;$i<scalar(@header);$i++){
	$heads{$header[$i]}=$i+1;
    }
    my @colnames=split /,/,$get_cols_by_header,-1;
    $col_ranges="";
    for my $cn (@colnames){
	if (exists $heads{$cn}){
	    $col_ranges.=$heads{$cn}.",";
	}
    }
    chop $col_ranges;
}
if ($skip_col>0){ $col_ranges="1-$skip_col,".$col_ranges }

@cols = &parseRanges($col_ranges);


my $data_line=0;
my @tuple;
my @tuple_all;
my $value;
my $beg;
my $end;
my $fixed_cols=0;
my $cols_ref;
my $max_cols = 0;
while($first_line ne "" or $_=<$file_ref>)
{
  if ($first_line ne ""){
    $_=$first_line;
    $first_line="";
  }

  $line++;

  if ($headers >= $line)
  {
      print;
  }
  else
    {
      if(defined($col_ranges))
	{
	  my $num_cols = &numTokens($delim_in);
	  if($num_cols > $max_cols or ($tight and $num_cols != $max_cols))
	    {
	      $max_cols = $num_cols;
	      @cols     = &parseRanges($col_ranges, $max_cols);
	      # print STDERR join(',',@cols), "\n";
	      if (!$zerobased){
		for(my $i = 0; $i < scalar(@cols); $i++)
		  {
		    $cols[$i]--;
		  }
	      }
	    }
	}
      
      if(/\S/ and not(/^\s*#/))
	{
	  if($multiple_delim_ins)
	    {
	      @tuple_all = split(/[$delim_in]+/,$_,-1);
	    }
	  else
	    {
	      @tuple_all = split($delim_in,$_,-1);
	    }
	  if($#tuple_all >= 0)
	    { chomp($tuple_all[$#tuple_all]); }
	  
	  @tuple=();
	  
	  if (not $invert)
	  {
	      foreach my $i (@cols)
	      {
		  if($preserve_empty_values or $i <= $#tuple_all)
		  {
		      push(@tuple,$tuple_all[$i]);
		  }
	      }
	  }
	  else
	  {
	      my %cols_hash;
	      for my $i (@cols)
	      {
		  $cols_hash{$i}=1;
	      }
	      for my $i (0..$#tuple_all)
	      {
		  if($i <= $#tuple_all and not $cols_hash{$i})
		  {
		      push(@tuple,$tuple_all[$i]);
		  }
	      }
	  }

	  my $result = join($delim_out,@tuple);

	  if(not($suppress_blanks) or $result =~ /\S/)
	    {
	      print $result, "\n";
	    }
	}
      elsif(not($suppress_blanks) or /\S/)
	{
	  print;
	}
    }
}

exit(0);

__DATA__

syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/cut.pl [OPTIONS] TAB_FILE

TAB_FILE is any tab-delim_inited file.  Can also be passed into standard
   input.

OPTIONS are:

-q: Quiet mode: turn verbosity off (default verbose)
-h NUM: Set the number of row headers to NUM (default 0).
-d DELIM: Change the input and output delimiter to DELIM (default <tab>).
-di DELIM: Change the input delimiter to DELIM (default <tab>).
-do DELIM: Change the output delimiter to DELIM (default <tab>).
-s:        Suppress blanks
-i:        Invert. Use the complement of the fields specified in -f.
-t:        Tight.  Tell the script to expect different number of columns in
           each row so that it needs to recompute the column boundaries for
           each row.  Note this option slows the script down somewhat.
-f RANGES: specify column ranges to include.  RANGES are comma-
           seperated lists of single columns or a range of columns
           for example:

                   5-6,2,1-3

           would select columns 1 through 6 except column 4.  Note
           that 2 is redundantly specified by no error results.
-n <str>:  specify columns by names (instead of by numbers with -f). Assumes
           the first line is a header, and translates column names in <str>
           into their respective column numbers in the header. Supports
           multiple ranges. Column names are assumed to be unique and not to
	   contain any commas.
-file <str>: get column names (as in -n) from specified file. each column name
           should appear on a separate line.
-0:        zero-based column numbers.
-e:        preserve empty values.
-sk:       number of header columns (default 0)

