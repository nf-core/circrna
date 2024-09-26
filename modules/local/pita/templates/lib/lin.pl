#!/usr/bin/perl

use strict;

my $arg='';
my $delim="\t";
my $i=1;
my $after=0;
my $neg = '';
my $inc=1;

while(@ARGV)
  {
    $arg = shift @ARGV;
    if($arg =~ /-([-]*\d+)/)
      {
        # Start counting from the number they pass in:
        $i=int($1);
      }
    elsif($arg eq '-d')
      {
        $delim = shift @ARGV;
      }
    elsif($arg eq '-a')
      {
        $after=1;
      }
    elsif($arg eq '-n')
      {
        $neg = '-';
      }
    elsif($arg eq '-i')
      {
        $inc = shift @ARGV;
      }
  }

# Add line numbers to the file.
while(<STDIN>)
  { 
     chop;
     if(not($after))
       { print STDOUT $neg, "$i", $delim, "$_\n"; }
     else
       { print STDOUT "$_", $delim, $neg, "$i\n"; }

     $i += $inc;
  }
