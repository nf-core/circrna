#!/usr/bin/perl

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libtable.pl";

use strict;

################# BEGIN MAIN ###########################

my $arg='';
my (@key1) = (1);
my (@key2) = (1);
my ($beg,$end) = (0,0);
my ($i,$j) = (0,0);
my ($file1,$file2) = ('','');
my ($key)='';
my (%values)={};
my (%exists)={};
my $delim = "\t";
my ($delim_in1, $delim_in2) = ("\t","\t");
my ($delim_out)  = ("\t");
my $delim_syn = "\t";
my (@tmp) = ();
my ($value1, $value2) = ('','');
my ($printable, $printable1, $printable2) = ('','','');
my ($suppress1, $suppress2, $suppressk) = (0,0,0);
my $negate=0;
my ($numeric) = 0;
my $empty = '';
my $fill = '';
my ($outer) = 0;
my ($outer_id) = 0;
my ($reverse) = 0;
my ($uppercase) = 0;
my $no_ladder=0;
my $hit = 0;
my $tmp;
my @fill_lines;
my $merge=0;
my $syn_file = "$ENV{JOIN_SYNONYMS}";
my %syn_pairs;
my $syn_pair;
my @syns_remaining;
my @syn_list;
my %syn_seen;
my $verbose=1;
my $header=0;
my $skip = 0;
my $ignore_missing_file=0;

while(@ARGV)
{
    $arg = shift @ARGV;
    if($arg eq '--help')
    {
      print STDOUT <DATA>;
      exit(0);
    }
    elsif($arg eq '-q')
    {
      $verbose = 0;
    }
    elsif($arg eq '-1')
      {
        $arg = shift @ARGV;
        @key1 = &parseRanges($arg);
      }
    elsif($arg eq '-2')
      {
        $arg = shift @ARGV;
        @key2 = &parseRanges($arg);
      }
    elsif($arg eq '-o' or $arg eq '--outer')
    {
      $fill = shift @ARGV;
      $outer = 1;
    }
    elsif($arg eq '-imf')
    {
      $ignore_missing_file = 1;
    }
    elsif($arg eq '-of')
    {
      $arg = shift @ARGV;
      if(not(open(FILE, $arg)))
      {
        print STDERR "Could not open file '$arg' to find outer text, skipping.\n";
      }
      else
      {
        while(<FILE>)
        {
          if(/\S/)
          {
            chomp;
            push(@fill_lines, $_);
          }
        }
      }
      $outer = 1;
    }
    elsif($arg eq '-oid')
    {
      $outer_id = 1;
    }
    elsif($arg eq '-skip')
    {
      $skip = shift @ARGV;
    }
    elsif($arg eq '-no_ladder')
    {
      $no_ladder = 1;
    }
    elsif($arg eq '-h1')
    {
      $header = 1;
    }
    elsif($arg eq '-h2')
    {
      $header = 2;
    }
    elsif($arg eq '-e' or $arg eq '--empty')
    {
      $empty = shift @ARGV;
    }
    elsif($arg eq '-m')
    {
      $merge = 1;
    }
    elsif($arg eq '-num')
      {
        $numeric = 1;
      }
    elsif($arg eq '-neg')
      {
        $negate = 1;
      }
    elsif($arg eq '-di1')
      {
        $delim_in1 = shift @ARGV;
      }
    elsif($arg eq '-di2')
      {
        $delim_in2 = shift @ARGV;
      }
    elsif($arg eq '-di')
      {
        $delim_in1 = shift @ARGV;
        $delim_in2 = $delim_in1;
      }
    elsif($arg eq '-do')
      {
        $delim_out = shift @ARGV;
      }
    elsif($arg eq '-ds')
    {
      $delim_syn = shift @ARGV;
    }
    elsif($arg eq '-syn')
    {
      $syn_file = shift @ARGV;
    }
    elsif($arg eq '-nosyn')
    {
      $syn_file = '';
    }

    # Suppress printing of values from table 1 (key will be printed however).
    elsif($arg eq '-s1')
    {
        $suppress1 = 1;
    }
    # Suppress printing of values from table 2 (key will be printed however).
    elsif($arg eq '-s2')
      {
        $suppress2 = 1;
      }
    elsif($arg eq '-sk')
      {
        $suppressk = 1;
      }
    elsif($arg eq '-r' or $arg eq '--reverse' or $arg eq '-rev')
      {
        $reverse = 1;
      }
    elsif($arg eq '-u')
      { $uppercase = 1; }

    elsif(length($file1)<1)
      {
        $file1 = $arg;
      }
    elsif(length($file2)<1)
      {
        $file2 = $arg;
      }
}

# open(SYN,">tmp_syns");

if(@fill_lines and $#fill_lines>=0)
{
  $fill = join($delim_out, @fill_lines);
}

if($reverse)
{
  $tmp = $suppress1;
  $suppress1 = $suppress2;
  $suppress2 = $tmp;
}

# See if the user supplied a synonyms file.  If so, extract synonyms for
# keys from the file and load it into a hash.
my %syns;
my $syn;
if(length($syn_file)>0)
{
  if(not(open(FILE, $syn_file)))
  {
    print STDERR "Could not open the synonyms file, '$syn_file' skipping.\n";
  }
  else
  {
    if($verbose)
      { print STDERR "Reading in synonyms from '$syn_file'..."; }
    while(<FILE>)
    {
      chomp;
      @tmp = split($delim_syn,$_,-1);
      if($uppercase)
      {
        for($i=0; $i<=$#tmp; $i++)
        {
          $tmp[$i] =~ tr/a-z/A-Z/;
        }
      }
      if($numeric)
      { 
        for($i=0; $i<=$#tmp; $i++)
        {
          $tmp[$i] =~ int($tmp[$i]);
        }
      }
      for($i=0; $i<=$#tmp; $i++)
      {
        if($tmp[$i] =~ /\S/)
        {
          for($j=0; $j<=$#tmp; $j++)
          {
            if($tmp[$j] =~ /\S/)
            {
              $syns{$tmp[$i]} .= $tmp[$j] . $delim_syn;
              $syn_pairs{$tmp[$i] . $delim_syn . $tmp[$j]} = 1;
            }
          }
        }
      }
    }
    # Post-process the synonyms: If there are synonyms of synonyms, then make
    # sure these are united etc..
    foreach $syn_pair (keys(%syn_pairs))
    {
      (@syns_remaining) = split($delim_syn,$syn_pair,-1);
      @syn_list=();
      # print SYN '[', join($delim_syn, @syns_remaining), "]: ";
      while(@syns_remaining)
      {
        $syn = shift @syns_remaining;
        if(not($syn_seen{$syn}))
        {
          $syn_seen{$syn} = 1;
          push(@syn_list,$syn);
          # print SYN "$syn ";
          @tmp = split($delim_syn, $syns{$syn},-1);
          # Add new synonyms to the list to be processed.
          for($j=0; $j<=$#tmp; $j++)
          {
            push(@syns_remaining, $tmp[$j]);
          }
        }
      }
      # print SYN "]\n";
      for($i=0; $i<=$#syn_list; $i++)
      {
        $syns{$syn_list[$i]} = '';
        for($j=0; $j<=$#syn_list; $j++)
        {
          $syns{$syn_list[$i]} .= $syn_list[$j];
          if($j<$#syn_list)
          {
            $syns{$syn_list[$i]} .= $delim_syn;
          }
        }
        # print SYN "[$syn_list[$i]] <=> [$syns{$syn_list[$i]}]\n";
      }
    }

    close(FILE);
    if($verbose)
      { print STDERR " done.\n"; }
  }
}

# $arg = $syns{'F27B3.1'};
# print "[$arg]\n";

if(length($file1)<1 or length($file2)<1)
{
  print STDERR <DATA>;
  exit(1);
}

if ($ignore_missing_file and !(-e $file1) and ($file1 ne "-")){
  open(F2,$file2);
  while(my $line=<F2>){
    my %check;
    chomp $line;
    my @tmp=split /$delim_in2/,$line,-1;
    my @out;
    for my $i (@key2){
      push @out,$tmp[$i-1];
      $check{$i-1}=1;
    }
    for (my $i=0;$i<scalar(@tmp);$i++){
      if(!$check{$i}){push @out,$tmp[$i]}
    }
    print join($delim_out,@out),"\n";
  }
  close(F2);
  exit;
}
if ($ignore_missing_file and !(-e $file2) and ($file2 ne "-")){
  open(F1,$file1);
  while(my $line=<F1>){
    my %check;
    chomp $line;
    my @tmp=split /$delim_in1/,$line,-1;
    my @out;
    for my $i (@key1){
      push @out,$tmp[$i-1];
      $check{$i-1}=1;
    }
    for (my $i=0;$i<scalar(@tmp);$i++){
      if(!$check{$i}){push @out, $tmp[$i]}
    }
    print join($delim_out,@out),"\n";
  }
  close(F1);
  exit;
}



for($i=0; $i<=$#key1; $i++)
  { $key1[$i]--; }

for($i=0; $i<=$#key2; $i++)
  { $key2[$i]--; }

@key1 = sort {$a <=> $b} @key1;
@key2 = sort {$a <=> $b} @key2;

# print STDERR "Key 1: [", join(',', @key1), "]\n",
#         "Key 2: [", join(',', @key2), "]\n",
#         "Input delimiter 1: [$delim_in1]\n",
#         "Input delimiter 2: [$delim_in2]\n",
#         "Output delimiter 1: [$delim_out1]\n",
#         "Output delimiter 2: [$delim_out2]\n",
#         "\n",
#         ;

# Read in the key-printable pairs from the second file:
my ($loops) = 0;
my ($passify) = 10000;
my $file_ref;
if($file2 =~ /\.gz$/)
{
  open(FILE, "zcat < $file2 |") or die("Could not open file '$file2'.");
  $file_ref = \*FILE;
}
elsif($file2 eq '-')
{
  $file_ref = \*STDIN;
}
else
{
  open(FILE, $file2) or die("Could not open file '$file2'.");
  $file_ref = \*FILE;
}
if($verbose)
{ 
  print STDERR "Reading relations from ", $file2 eq '-' ? "standard input" 
                  : "file $file2";
}
my $header_data='';
my $line=0;
while(<$file_ref>)
{
    if(/\S/ and (!(/^\s*#/) or $no_ladder))
      {
	$line++;
	if($line==1 and $header==2)
	  { $header_data = $_; }
	chomp;
        @tmp = split($delim_in2,$_,-1);
# print STDERR "\n2: tmp: [", join('|',@tmp), "]\n";
# print STDERR "2: key cols: [", join('|',@key2), "]\n";
        $key='';
        for($i=$#key2; $i>=0; $i--)
          {
            my $key_part = splice(@tmp,$key2[$i], 1);
	    # $key .= length($key)>0 ? ($delim_out . $key_part) : $key_part;
	    $key = length($key)>0 ? ($key_part . $delim_out . $key) : $key_part;
# print STDERR "1: key: [$key]\n";
# print STDERR "Before splice: [", join('|',@tmp), "] [$key]\n";
# print STDERR "After splice: [", join('|',@tmp), "] [$key]\n";
            # $key .= splice(@tmp, $i-1, 1) . $delim_out;
          }
        # Get rid of the last delimiter we added:
        if($numeric)
          { $key = int($key); }
        if($uppercase)
          { $key =~ tr/a-z/A-Z/; }
# print STDERR "2: tmp: [", join('|',@tmp), "]\n";
# print STDERR "2: key: [$key]\n";
        $tmp = join($delim_out, @tmp);
	if ($#tmp<0){$tmp="--BLANKNOPRINT--"}
        $values{$key} = $tmp;
        $exists{$key} = 1;

        # Record this value with all the synonym keys as well.

        @tmp = split($delim_syn, $syns{$key},-1);
        # foreach $key (@tmp)
        for($j=0; $j<=$#tmp; $j++)
        {
          if($numeric)
            { $key = int($key); }
          if($uppercase)
            { $key =~ tr/a-z/A-Z/; }
          $values{$key} = $tmp;
          $exists{$key} = 1;
        }

        $loops++;
        if($verbose and $loops%$passify==$passify-1)
          {
            print STDERR '.';
          }
      }
}
if($verbose)
  { print STDERR " done.\n"; }
close($file_ref);

# $arg = $syn{'F15G10.1'};
# print STDERR "[$arg]\n";

# Read in the key-printable pairs from the first file and print out
# the joined key:
my $file_ref;
if($file1 =~ /\.gz$/)
{
  open(FILE, "zcat < $file1 |") or die("Could not open file '$file1'.");
  $file_ref = \*FILE;
}
elsif($file1 eq '-')
{
  $file_ref = \*STDIN;
}
else
{
  open(FILE, $file1) or die("Could not open file '$file1'.");
  $file_ref = \*FILE;
}
if($verbose)
  { print STDERR "Joining on file $file1"; }
$loops=0;
my $found;
while(<$file_ref>)
  {
    if ($skip > $loops)
      {
	$loops++;
	print;
      }
    elsif(/\S/ and (!(/^\s*#/) or $no_ladder))
      {
	chomp;
        @tmp = split($delim_in1,$_,-1);
        $key='';
        for($i=$#key1; $i>=0; $i--)
          {
            my $key_part = splice(@tmp,$key1[$i], 1);
	    $key = length($key)>0 ? ($key_part . $delim_out . $key) : $key_part;
          }
        # Get rid of the last delimiter we added:
        if($numeric)
          { $key = int($key); }
        if($uppercase)
          { $key =~ tr/a-z/A-Z/; }

        $value1 = join($delim_out, @tmp);
	if ($#tmp<0){$value1="--BLANKNOPRINT--"}
        $value2 = $values{$key};
        $found = $exists{$key};

# print STDERR "1: key: [$key]\n";
# print STDERR "1: value1: [$value1]\n";
# print STDERR "1: value2: [$value2]\n";
# print STDERR "1: found: [$found]\n";

        # See if this key matches any of the key's synonyms
        if(not($found))
        {

          @tmp = split($delim_syn, $syns{$key},-1);
          while(@tmp and not($found))
          {
            $tmp = shift @tmp;
            $found = $exists{$tmp};
            if($found)
            {
              $value2 = $values{$tmp};
            }
          }
        }

        if((not($negate) and $found))
          {
            $hit = 1;
         }
	elsif($negate and not($found))
	  {
	    $hit=1 ;
	    $value2 = "--BLANKNOPRINT--";
	  }
        elsif($outer)
          {
            $value2 = $fill;
	    if ($value2 eq ""){$value2 = "--BLANKNOPRINT--"}
            $hit = 1;
          }
        elsif($outer_id)
	  {
            $value2 = $key;
            $hit = 1;
	  }
        else
          {
            $hit = 0;
          }

        if($merge & $hit)
        {
          $value1 = $value2;
          $value2 = "--BLANKNOPRINT--";
        }
        elsif($merge & not($hit))
        {
          $hit = 1;
          $value2 = "--BLANKNOPRINT--";
        }

        if($hit)
          {
            # Swap the two values if we're supposed to print the second value
            # before the first value.
            if($reverse)
            {
              $tmp = $value1;
              $value1 = $value2;
              $value2 = $tmp;
            }

	    my @printed;

	    if(!$suppressk){push @printed,$key}
	    if(!$suppress1 and $value1 ne "--BLANKNOPRINT--"){push @printed,length($value1)<1?$empty:$value1}
	    if(!$suppress2 and $value2 ne "--BLANKNOPRINT--"){push @printed,length($value2)<1?$empty:$value2}
	    print join($delim,@printed),"\n";
          }

        $loops++;
        if($verbose and $loops%$passify==$passify-1)
          {
            print STDERR '.';
          }
      }
  }
if($verbose)
  { print STDERR " done.\n"; }
close(FILE);

exit(0);

################# END MAIN #############################

__DATA__
syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/join.pl [OPTIONS] FILE1 FILE2

This script takes two tables, contained in delimited files, as input and
produces a new table that is a join of FILE1 and FILE2.  The script assumes
the keys can be found in the first column of the two files and that columns
are delimited by single tabs (both of these assumptions can be changed, see
OPTIONS).  If FILE1 contains the tuple (A, W, X) and FILE2 contains the
tuple (A, Y, Z) then (A, W, X, Y, Z) will be contained in the output.
NOTE: FILE2 will be loaded into memory.

OPTIONS are:

-m: Merge - if key exists in file 2, use value from file 2 else
    use the value from file 1.
-q: Quiet mode: turn verbosity off (default is verbose)

-1 COL: Include column COL from FILE1 as part of the key (default is 1).
         Multiple columns may be specified in which case keys are constructed
         by concaternating each column in the order specified (delimiting
         character used equal to the output delimiter, see -do flag).

-2 COL: Include column COL from FILE2 as part of the key (default is 1).
         See -1 option for discussion of multiple columns.

-skip NUM: Skip the first NUM lines of FILE1 and print them as they are (first line
            is 1, default value is 0).
-o FILL: Do a left outer join.  If a key in FILE1 is not in FILE2, then the
--outer          tuple from FILE1 is printed along with the text in FILL in place of
          a tuple from FILE2 (by default these tuples are not reported in the
          result).  See -of option also to supply FILL from a file, and -oid to replace
          FILL with the id from FILE1.

-of FILE: Same as -o, but use the text in FILE as the text to FILL.

-oid: Same as -o, but use the id from FILE1 as the text to FILL.

-e EMPTY: Set the empty string to EMPTY (default is blank).  If both keys
           exist in FILE1 and FILE2 but one tuple is blank, then the empty
           character EMPTY will be printed.

-num: Treat the keys as numeric quantities (default is off: treat keys as
        strings).  If this is turned on, each key will be forced into an
        integer quantity.

-neg: Negative output -- print out the keys from FILE1 that are not in FILE2.
        These are equivalent to those keys that would be left out of the join
        or those that would have a FILL tuple in a left outer join (see -o
        option).

-di1 DELIM: Set the input delimiter for FILE1 to DELIM (default is tab).

-di2 DELIM: Set the input delimiter for FILE2 to DELIM (default is tab).

-di DELIM: Set the input delimiters for both FILE1 and FILE2 to DELIM (default
        is tab).  Equivalent to using both the -di1 and -di2 options.

-do DELIM: Set the output delimiter to DELIM.  Note there is only one output
        delimiter (not two: one for FILE1 and FILE2); this forces the
        output to a common delimitation (default is tab).

-ds DELIM: Set the delimiter for synonyms to DELIM.  This is used for reading
        synonyms for keys (see the -syn option) (default is tab).

-syn FILE: Specify that synonyms for keys can be found in the file FILE.
        Each line in FILE should contain synonyms for only one key, seperated
        by the delimiter specified with the -ds option).

-nosyn: Ignore synonym file in the JOIN_SYNONYMS environment variable if it
        exists.

-s1: Suppress printing of tuples from FILE1.  The key is printed, followed by
        the tuple found in FILE2.

-s2: Suppress printing of tuples from FILE2.

-sk: Suppress printing of keys.

-r: Reverse -- Instead of printing <key> <tuple1> <tuple2>, print
        <key> <tuple2> <tuple1> where <key> is the key, <tuple1> is a tuple
        found from FILE1 and <tuple2> is a tuple found in FILE2.

-u: Uppercase.  For non-numeric keys, force any letter to be uppercase before
        attempting to do key lookups.  Keys from both FILE1 and FILE2 will be
        converted to uppercase before attempting the join.

-imf: Ignore missing file. If file1 does not exist, simply print file2 (and
        vice versa).
-no_ladder: Turn off default behavior of ignoring lines with #

