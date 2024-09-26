#!/usr/bin/perl

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";

use strict;

my $beg          = undef;
my $end          = undef;
my $count_blanks = 1;
my $quiet        = 0;
my @skip_lines;
my @select_lines;
my $fin = \*STDIN;
my $fraction=0;
while(@ARGV)
{
  my $arg = shift @ARGV;
  if($arg eq '--help')
  {
    print STDOUT <DATA>;
    exit(0);
  }
  elsif($arg eq '-f')
  {
    $fraction = 1;
  }
  elsif($arg eq '-b')
  {
    $count_blanks = 0;
  }
  elsif(not(defined($beg)))
  {
    $beg = $fraction?$arg:int($arg);
  }
  elsif(not(defined($end)))
  {
    $end = $fraction?$arg:int($arg);
  }
  elsif($arg eq '-skip')
  {
     my $lines_str = shift @ARGV;
     @skip_lines = sort {$a <=> $b} parseRanges($lines_str);
  }
  elsif($arg eq '-select')
  {
     my $lines_str = shift @ARGV;
     @select_lines = sort {$a <=> $b} parseRanges($lines_str);
  }
  elsif($arg eq '-quiet')
  {
     $quiet = 1;
  }
  elsif(-f $arg or -l $arg)
  {
     open($fin, $arg) or die("Could not open file '$arg'");
  }
  else
  {
    die("/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/body.pl: Bad argument '$arg' given.  Use --help for help.");
  }
}

if(not(defined($beg)))
{
  $beg = 1;
}

if(not(defined($end)))
{
  $end = -1;
}

my $num_lines = undef;
my $tmp_file  = undef;
if($end < -1 or ($beg<1 and $beg>0) or ($end>-1 and $end<0) or ($end<1 and $end>0))
{
   $tmp_file = 'tmp_' . time . '.' . rand() . './nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/body.pl';
   open(TMP, ">$tmp_file") or die("Could not open temporary file '$tmp_file' for writing");
   while(<$fin>)
     { print TMP; }
   close(TMP);

   my $wc = `wc -l $tmp_file`;
   my @tuple = split(/\s+/,$wc);
   $num_lines = $tuple[0];
   print STDERR "The file has $num_lines number of lines.\n" unless ( $quiet );

   open($fin, "<$tmp_file") or die("Could not open the temporary file '$tmp_file' for reading");

}

my $line = 0;
my $skip_counter = 0;
my $select_counter = 0;
if(defined($num_lines)){
  if ($end<-1){
    $end=$num_lines + $end + 1;
  }
  if (($end>-1 and $end<0)){
    $end = $num_lines + sprintf("%.f",$end*$num_lines);
  }
  if (($end<1 and $end>0)){
    $end = sprintf("%.f",$end*$num_lines);
  }
  if ($beg<1 and $beg>0){
    $beg = sprintf("%.f",$beg*$num_lines) + 1;
  }
}
while(<$fin>)
{
  if($count_blanks or /\S/)
  {
    $line++;
    if (@skip_lines >= $skip_counter and $line == $skip_lines[$skip_counter])
    {
       $skip_counter++;
       next;
    }
    if (@select_lines>0 and $line != $select_lines[$select_counter])
    {
      next;
    }
    $select_counter++;

 
    if($line >= $beg and ($end == -1 or $line <= $end))
    {
      print;
    }
  }
}

if(defined($tmp_file))
{
   system("rm -f $tmp_file");
}

exit(0);

__DATA__
Syntax: /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/body.pl BEG END < FILE

BEG, END are the beginning and end lines (inclusive) to select from
the file.  If END=-1 then the rest of the file is included for example BEG=2 END=-1
returns the whole file except the first row.

OPTIONS are:

-b:               Do *not* include blank lines when counting (default counts them).
-skip <n1,n2...>: Exclude line numbers n1,n2...
-select <n1,n2>:  Select line numbers n1,n2...
-quiet:           Do not print any message to STDERR.
-f:               Allow giving BEG and END as fractions of number of lines in file (for example /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/body.pl 4 -0.2)


