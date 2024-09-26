#!/usr/bin/perl

use strict;

my $arg;
my $verbose=1;

while(@ARGV)
{
  $arg = shift @ARGV;

  if($arg eq '--help')
  {
    &printSyntax();
    exit(0);
  }
  elsif($arg eq '-q')
  {
    $verbose = 0;
  }
  else
  {
    &printSyntax();
    exit(1);
  }
}

if(not(-T STDIN))
{
  if($verbose)
  {
    print STDERR "File is not a text file, skipping.\n";
  }
  exit(1);
}

while(<STDIN>)
{
  s///g;
  print;
}

exit(0);

sub printSyntax
{
  print STDERR "Syntax: $0 [OPTIONS] < FILE\n",
	"\n",
	"where FILE is any text file",
	"\n",
	"where OPTIONS are\n",
	"\n",
	"-q: Quiet mode (default is verbose)\n",
  	"\n";
}

