#!/usr/bin/perl

use strict;

sub format_number
{
  my $number = $_[0];
  my $sig_digits = $_[1];

  my @num = split(/\./, $number);

  my $suffix = "";
  if ($num[1] =~ /(.*)e(.*)/) { $num[1] = $1; $suffix = "e$2"; }

  my $decimal = "";

  if ($num[0] == 0)
  {
    my $lsb_len = length($num[1]);

    $num[1] =~ /([1-9].*)/;

    my $significant_len = length($1);

    $decimal = substr($num[1], 0, $lsb_len - $significant_len + $sig_digits);
  }
  else
  {
    $decimal = substr($num[1], 0, $sig_digits);
  }

  if (length("$decimal$suffix") > 0) { return "$num[0].$decimal$suffix"; }
  else { return "$num[0]"; }
}

#print format_number($ARGV[0],$ARGV[1]) . "\n";

1
