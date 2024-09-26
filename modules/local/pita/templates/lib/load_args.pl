#!/usr/bin/perl

#---------------------------------------------------------
# load_args
#---------------------------------------------------------
sub load_args (\@)
{
  my @args = @{@_[0]};

  my %result;

  my $num_args = @args;

  for (my $current_arg = 0; $current_arg < $num_args; $current_arg++)
  {
    my $arg = $args[$current_arg];

    if ($arg =~ /-([^\s]+)/)
    {
      my $arg_name = $1;

      my $next_arg = $args[$current_arg + 1];

      if ($next_arg =~ /^-([^\s])/ || length($next_arg) == 0)
      {
	$result{$arg_name} = "1";
      }
      else
      {
	$next_arg =~ s/[\']//g;
	$next_arg =~ s/[\"]//g;

	$result{$arg_name} = $next_arg;
	$current_arg++;
      }

      #print "load_args: result{$arg_name}=$result{$arg_name}\n";
    }
  }

  return %result;
}

#---------------------------------------------------------
# get_arg
#---------------------------------------------------------
sub get_arg ($$\%)
{
  my ($arg, $default, $str_args) = @_;
  my %args = %$str_args;

  if (length($args{$arg}) > 0) { return $args{$arg}; }
  else { return $default; }
}

#---------------------------------------------------------
# echo_arg
#---------------------------------------------------------
sub echo_arg ($$\%)
{
  my ($arg, $str_args) = @_;
  my %args = %$str_args;

  if (length($args{$arg}) > 0) { return "-$arg $args{$arg} "; }
  else { return ""; }
}

sub echo_arg_quoted ($$\%)
{
  my ($arg, $str_args) = @_;
  my %args = %$str_args;

  if (length($args{$arg}) > 0) { return "-$arg \"$args{$arg}\" "; }
  else { return ""; }
}
#---------------------------------------------------------
# echo_arg_equals
#---------------------------------------------------------
sub echo_arg_equals ($$\%)
{
  my ($arg, $str_args) = @_;
  my %args = %$str_args;

  if (length($args{$arg}) > 0) { return "$arg" . "=" ."\"$args{$arg}\" "; }
  else { return ""; }
}
#---------------------------------------------------------
# get_extended_arg
#---------------------------------------------------------
sub get_extended_arg ($\%)
{
  my ($arg, $str_args) = @_;
  my %args = %$str_args;

  my @res;
  my $id = 1;
  my $done = 0;
  while ($done == 0)
  {
      my $extended_arg = &get_arg("$arg$id", "", \%args);
      if (length($extended_arg) > 0)
      {
	  push(@res, $extended_arg);

	  $id++;
      }
      else
      {
	  $done = 1;
      }
  }

  return @res;
}

#---------------------------------------------------------
# get_full_arg_command
#---------------------------------------------------------
sub get_full_arg_command (\@)
{
  my @args = @{@_[0]};

  my $res = "";

  for (my $current_arg = 0; $current_arg < @args; $current_arg++)
  {
    my $arg = $args[$current_arg];

    if ($arg =~ /-([^\s]+)/)
    {
      $res .= "$arg ";
    }
    else
    {
    	$res .= "'$arg' ";
    }
  }

  return $res;
}

1
