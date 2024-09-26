use strict;

#------------------------------------------------------------------------
# $int getNumCols($string file, $string delim="\t")
#------------------------------------------------------------------------
sub getNumCols
{
   my ($file, $delim) = @_;
   $delim = not(defined($delim)) ? "\t" : $delim;
   my $fp;
   open($fp, $file) or die("Could not open file '$file'");
   my $done = 0;
   my $num_cols = undef;
   while(not(defined($num_cols)) and not(eof($fp)))
   {
      my $line  = <$fp>;
      if($line =~ /\S/)
      {
         my @tuple = split($delim, $line);
         $num_cols = scalar(@tuple);
         # print STDERR "$num_cols : ", join("-", @tuple), "'\n";
      }
   }
   return $num_cols;
}

#------------------------------------------------------------------------
# ($string, \@list) splitKeyAndValue ($string line, \@list key_cols,
#                                     \@list sorted_key_cols, $string delim)
#------------------------------------------------------------------------
sub splitKeyAndValue
{
  my ($line,$key_cols,$sorted_key_cols,$delim) = @_;
  $sorted_key_cols = defined($sorted_key_cols) ? $sorted_key_cols : $key_cols;
  $delim = defined($delim) ? $delim : "\t";
  my @tuple = split($delim, $line);
  my @tuple_copy = @tuple;
  
  if($#tuple == -1)
    { return (undef,undef); }

  chomp($tuple[$#tuple]);
  chomp($tuple_copy[$#tuple]);

  # Remove key entries from the tuple in reverse order to
  # preserve the indexing.
  my @key;
  my $j=0;
  for(my $i=scalar(@{$sorted_key_cols})-1; $i>=0; $i--)
  {
    $j++;
    push(@key, $tuple_copy[$$key_cols[$j]]);
    splice(@tuple,$$sorted_key_cols[$i], 1);
  }
  return (join($delim,@key),\@tuple);
}

#---------------------------------------------------------------------------
# \@list listRead ($string file, $string delim="\t", \@list key_cols=(1))
#---------------------------------------------------------------------------
sub tableRead
{
   my ($file, $delim, $key_cols) = @_;
   $delim = not(defined($delim)) ? "\t" : $delim;
   if(not(defined($key_cols)))
   {
      my @key_cols = (1);
      $key_cols = \@key_cols;
   }
   my @sorted_key_cols = sort { $a <=> $b; } @{$key_cols};
   my $filep;
   open($filep, $file) or die("Could not open file '$file'");
   my @tuples;
   while(<$filep>)
   {
      chomp;
      my $line = $_;
      my ($key, $tuple) = &splitKeyAndValue($line, $key_cols, \@sorted_key_cols, $delim);
      my @key_tuple = ($key);
      push(@key_tuple, @{$tuple});
      push(@tuples, \@key_tuple);
   }
   close($filep);

   return \@tuples;
}

#------------------------------------------------------------------------
# \@list emptyTuple ($int n)
#------------------------------------------------------------------------
sub emptyTuple
{
  my $n = shift;
  my @blanks = ();

  for(my $i=0; $i<$n; $i++)
    { push(@blanks,''); }
  
  return \@blanks;
}

#------------------------------------------------------------------------
# $int getTupleLength ($string line, $string delim)
#------------------------------------------------------------------------
sub getTupleLength
{
   my ($line, $delim) = @_;
   $line  = defined($line)  ? $line  : $_;
   $delim = defined($delim) ? $delim : "\t";

   my @tuple  = split($delim, $line);
   return ($#tuple + 1);
}

#------------------------------------------------------------------------
# $int compareKeys ($string key1, $string key2)
#------------------------------------------------------------------------
sub compareKeys
{
  my ($key1, $key2) = @_;

  if(not(defined($key1)) or length($key1)==0)
    { return 1; }
  
  if(not(defined($key2)) or length($key2)==0)
    { return -1; }
  
  return($key1 cmp $key2);
}

#------------------------------------------------------------------------
# $string getSubTupleString ($string line, \@list cols, $string delim)
#------------------------------------------------------------------------
sub getSubTupleString
{
  my ($line,$cols,$delim) = @_;
  $delim = defined($delim) ? $delim : "\t";
  my @tuple = split($delim, $line);

  if($#tuple==-1 or scalar(@{$cols})<=0)
    { return undef; }

  chomp($tuple[$#tuple]);

  my $sub = $tuple[$$cols[0]];
  for(my $i=1; $i<scalar(@{$cols}); $i++)
    { $sub .= $delim . $tuple[$$cols[$i]]; } 
  return $sub;
}

#------------------------------------------------------------------------
# \@list getSubTuple ($string line, \@list cols, $string delim)
#------------------------------------------------------------------------
sub getSubTuple
{
  my ($line,$cols,$delim) = @_;
  $delim = defined($delim) ? $delim : "\t";
  my @tuple = split($delim, $line);
  my @sub_tuple;

  chomp($tuple[$#tuple]);

  for(my $i=0; $i<scalar(@{$cols}); $i++)
  {
    push(@sub_tuple,$tuple[$$cols[$i]]);
  }
  return \@sub_tuple;
}

#------------------------------------------------------------------------
# $string getTupleEntry ($string line, int $i, $string delim)
#------------------------------------------------------------------------
sub getTupleEntry # ($line,$i,$delim)
{
  my ($line,$i,$delim) = @_;
  my @cols = ($i);
  return &getSubTupleString($line,\@cols,$delim);
}

#------------------------------------------------------------------------
# void forceTuplePrecision (\@list tuple, $double precision, $int skip)
#------------------------------------------------------------------------
sub forceTuplePrecision
{
  my ($tuple,$precision,$skip) = @_;
  $skip = (not(defined($skip)) or $skip < 0) ? 0 : $skip;
  $precision = 10**$precision;
  for(my $i=$skip; $i<scalar(@{$tuple}); $i++)
  {
    if($$tuple[$i] =~ /\d/)
    {
      $$tuple[$i] = int($$tuple[$i]*$precision) / $precision;
    }
  }
}

#------------------------------------------------------------------------
# $string forcePrecision ($string value, $double precision)
#------------------------------------------------------------------------
sub forcePrecision
{
   my ($value, $precision) = @_;
   if($value =~ /\d/)
   {
     $value = int($value*(10**$precision)) / (10**$precision);
   }
   return $value;
}

#------------------------------------------------------------------------
# $int getTuplePrecision (\@list tuple, $int skip)
#------------------------------------------------------------------------
sub getTuplePrecision
{
  my ($tuple,$skip) = @_;
  $skip = (not(defined($skip)) or $skip < 0) ? 0 : $skip;
  my $precision = 0;
  for(my $i=$skip; $i<scalar(@{$tuple}); $i++)
  {
    my $num = ($$tuple[$i] =~ s/(\d)/\1/g);
    if($num > $precision)
      { $precision = $num; }
  }
  return $precision;
}

#------------------------------------------------------------------------
# (\@list lines, $string last_line)
# getLinesWithIdenticalKeys ($fp file, \@list key_cols,
#                            $string last_line, string $delim)
#------------------------------------------------------------------------
sub getLinesWithIdenticalKeys
{
  my ($file,$key_cols,$last_line,$delim) = @_;
  my $key       = undef;
  my $first_key = undef;
  my @lines     = ();
  my $done      = 0;

  if(eof($file))
  {
    @lines = ($last_line);
    return (\@lines,undef);
  }

  while(not($done))
  {
    my $line   = defined($last_line) ? $last_line : <$file>;
    $last_line = undef;
    $key       = &getSubTupleString($line,$key_cols,$delim);
    $first_key = defined($first_key) ? $first_key : $key;
    my $order  = &compareKeys($key,$first_key);
    if($order == 0)
      { push(@lines, $line); } 
    else
    { 
      $done = 1;
      $last_line = $line;
    }
    
    if(eof($file))
      { $done = 1; }
  }
  return (\@lines,$last_line);
}

##------------------------------------------------------------------------
## void sequentialJoin ($string out_file, @list in_files)
##------------------------------------------------------------------------
sub sequentialJoin
{
   my ($out_file, @in_files) = @_;
   my $in_file = $in_files[0];
   my $join = "cat $in_file ";
   for(my $i = 1; $i <= $#in_files; $i++)
   {
      my $in_file   = $in_files[$i];
      $join .= "| /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/join.pl -q - $in_file ";
   }
   $join .= "> $out_file";

   print STDERR "--> $join <--\n";
   `$join`;
   print STDERR "--> Done. <--\n";
}


1
