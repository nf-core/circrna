#!/usr/bin/perl

use strict;

my $MUTUAL_INFORMATION = "mi";
my $DOT_PRODUCT = "dot";
my $CORRELATION_COEFFICIENT = "cor";
my $LOG2 = log(2);

#--------------------------------------------------------------------------------
# DEBUG_VO
#--------------------------------------------------------------------------------
sub DEBUG_VO
{
#  print $_[0];
}

#---------------------------------------------------------
# vec_avg
#---------------------------------------------------------
sub vec_avg($)
{
  my $vec = shift;

  my $sum = vec_sum($vec);

  my $entries = vec_count_full_entries($vec);
  #print "sum=$sum entries=$entries\n";
  return $entries > 0 ? $sum / $entries : "";
}

#---------------------------------------------------------
# best_dot_product
#---------------------------------------------------------
sub vec_best_dot_product($$)
{
  my ($vecX, $vecY) = @_;

  my $best_dot = "";

  my $sizeX = @$vecX;
  my $sizeY = @$vecY;

  if ($sizeX >= $sizeY)
  {
      for (my $X_offset = 0; $sizeX - $X_offset >= $sizeY; $X_offset++)
      {
	  my $dot = 0;
	  for (my $i = 0; $i < $sizeY; $i++)
	  {
	      $dot += ($$vecX[$i + $X_offset] * $$vecY[$i]);
	  }
	  
	  if (length($best_dot) == 0 or $dot > $best_dot)
	  {
	      $best_dot = $dot;
	  }
      }
  }
  else
  {
      for (my $Y_offset = 0; $sizeY - $Y_offset >= $sizeX; $Y_offset++)
      {
	  my $dot = 0;
	  for (my $i = 0; $i < $sizeX; $i++)
	  {
	      $dot += ($$vecY[$i + $Y_offset] * $$vecX[$i]);
	  }
	  
	  if (length($best_dot) == 0 or $dot > $best_dot)
	  {
	      $best_dot = $dot;
	  }
      }
  }

  return $best_dot;
}

#---------------------------------------------------------
# center
#---------------------------------------------------------
sub vec_center($)
{
  my $vec = shift;

  my $avg = vec_avg($vec);
  my $std = vec_std($vec);

  my @res;
  for (my $i = 0; $i < @$vec; $i++)
  {
    $res[$i] = (length($$vec[$i]) > 0 and $std > 0) ? (($$vec[$i] - $avg) / $std) : "";
  }

  return @res;
}

#---------------------------------------------------------
# center (by reference)
#---------------------------------------------------------
sub vec_center_by_ref($)
{
  my $vec = shift;

  my $avg = vec_avg($vec);
  my $std = vec_std($vec);

  for (my $i = 0; $i < @$vec; $i++)
  {
    $$vec[$i] = (length($$vec[$i]) > 0 and $std > 0) ? (($$vec[$i] - $avg) / $std) : "";
  }

}

#---------------------------------------------------------
# vec_compute_score
#---------------------------------------------------------
sub vec_compute_score ($$$)
{
  my ($vecX, $vecY, $op) = @_;

  if ($op eq $MUTUAL_INFORMATION)         { return vec_mutual_information($vecX, $vecY); }
  elsif ($op eq $DOT_PRODUCT)             { return vec_dot_product($vecX, $vecY); }
  elsif ($op eq $CORRELATION_COEFFICIENT) { return vec_correlation($vecX, $vecY); }
}

#---------------------------------------------------------
# vec_correlation
#---------------------------------------------------------
sub vec_correlation ($$)
{
  my ($vecX, $vecY) = @_;

  my $total = @$vecX;

  my $sum_X = 0;
  my $sum_XX = 0;
  my $sum_Y = 0;
  my $sum_YY = 0;
  my $sum_XY = 0;
  my $num = 0;

  for (my $i = 0; $i < $total; $i++)
  {
         my $p_X = $$vecX[$i];
         my $p_Y = $$vecY[$i];

         if($p_X =~ /\S/ and $p_Y =~ /\S/)
         {
            $sum_X  += $p_X;
            $sum_XX += $p_X * $p_X;
            $sum_Y  += $p_Y;
            $sum_YY += $p_Y * $p_Y;
            $sum_XY += $p_X * $p_Y;
            $num++;
         }
  }

  my $correlation = undef;

  my $numerator = ($num * $sum_XY) - ($sum_X * $sum_Y);
  my $var_X = ($num * $sum_XX) - ($sum_X * $sum_X);
  my $var_Y = ($num * $sum_YY) - ($sum_Y * $sum_Y);

  if ($var_X > 0 && $var_Y > 0)
  {
      my $denominator = sqrt($var_X * $var_Y);
      $correlation = $numerator / $denominator;
  }

  return $correlation;
}

#---------------------------------------------------------
# vec_count_full_entries
#---------------------------------------------------------
sub vec_count_full_entries($)
{
  my $vec = shift;

  my $entries = 0;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0) { $entries++; } }

  return $entries;
}

#---------------------------------------------------------
# dot_product
#---------------------------------------------------------
sub vec_dot_product ($$)
{
  my ($vecX, $vecY) = @_;

  my $dot = 0;
  for (my $i = 0; $i < @$vecX; $i++)
  {
      $dot += ($$vecX[$i] * $$vecY[$i]);
  }
  return $dot;
}

#---------------------------------------------------------
# dot_product_normedX
#---------------------------------------------------------
sub vec_dot_product_normedX ($$)
{
  my ($vecX, $vecY) = @_;

  my $dot = 0;
  my $norm = 0;
  for (my $i = 0; $i < @$vecX; $i++)
  {
      $dot += ($$vecX[$i] * $$vecY[$i]);
      $norm += abs($$vecX[$i]);
  }

  return ($norm > 0) ? ($dot/$norm) : "";
}

#---------------------------------------------------------
# entropy
#---------------------------------------------------------
sub vec_entropy ($)
{
  my $vecX = shift;

  DEBUG_VO("entropy for\n");
  DEBUG_VO("   Vec: @$vecX\n");

  my $total = @$vecX;

  my %COUNTS_X;
  for (my $i = 0; $i < $total; $i++)
  {
      my $key_X = $$vecX[$i];

      if (length($COUNTS_X{$key_X} == 0)) { $COUNTS_X{$key_X} = 1; } else { $COUNTS_X{$key_X} = $COUNTS_X{$key_X} + 1; }
  }

  my $entropy = 0;
  foreach my $key_X (keys %COUNTS_X)
  {
      my $p_X = $COUNTS_X{$key_X} / $total;

      DEBUG_VO("      PX($key_X) = $p_X\n");

      $entropy -= ($p_X) * log($p_X) / $LOG2;
  }

  return $entropy;
}

#---------------------------------------------------------
# intersect
#---------------------------------------------------------
sub vec_intersect ($$)
{
  my ($vecX, $vecY) = @_;

  my @res;
  my $counter = 0;

  my %h1;
  for (my $i = 0; $i < @$vecX; $i++) { $h1{$$vecX[$i]} = "1"; }

  for (my $i = 0; $i < @$vecY; $i++) { if ($h1{$$vecY[$i]} eq "1") { $res[$counter++] = $$vecY[$i]; $h1{$$vecY[$i]} = ""; } }

  return @res;
}

#---------------------------------------------------------
# kl distance
#---------------------------------------------------------
sub vec_kl_distance ($$$)
{
  my ($vecX, $vecY, $variable_dimension) = @_;

  my $res = 0;

  for (my $i = 0; $i < @$vecX; $i += $variable_dimension)
  {
      for (my $j = 0; $j < $variable_dimension; $j++)
      {
	  $res += $$vecX[$i + $j] * log($$vecX[$i + $j] / $$vecY[$i + $j]);
	  $res += $$vecY[$i + $j] * log($$vecY[$i + $j] / $$vecX[$i + $j]);
      }
  }

  return $res;
}

#---------------------------------------------------------
# vec_min
#---------------------------------------------------------
sub vec_min ($)
{
  my $vec = shift;

  my $min = 1e300;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] < $min) { $min = $$vec[$i]; } }

  return $min;
}

#---------------------------------------------------------
# vec_first_entry_ge_ind
#---------------------------------------------------------
sub vec_first_entry_ge_ind($)
{
  my $vec = shift;
  my $val = $_[1];

  my $ind = "";
  for (my $i = 0; $i < @$vec; $i++) {
    if ($$vec[$i]>=$val){$ind=$i;last;}
  }

  return $ind;
}

#---------------------------------------------------------
# vec_first_entry_le_ind
#---------------------------------------------------------
sub vec_first_entry_le_ind($)
{
  my $vec = shift;
  my $val = $_[1];

  my $ind = "";
  for (my $i = 0; $i < @$vec; $i++) {
    if ($$vec[$i]<=$val){$ind=$i;last;}
  }

  return $ind;
}


#---------------------------------------------------------
# vec_first_nonempty_entry_ind
#---------------------------------------------------------
sub vec_first_nonempty_entry_ind ($)
{
  my $vec = shift;

  my $ind = "";
  for (my $i = 0; $i < @$vec; $i++) {
    if ($$vec[$i] ne ""){$ind=$i;last;}
  }

  return $ind;
}

#---------------------------------------------------------
# vec_max_ind
#---------------------------------------------------------
sub vec_max_ind ($)
{
  my $vec = shift;
  my $ind = "";
  my $max = -1e300;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] > $max) { $max = $$vec[$i]; $ind = $i} }

  return $ind;
}

#---------------------------------------------------------
# vec_min_ind
#---------------------------------------------------------
sub vec_min_ind ($)
{
  my $vec = shift;
  my $ind = "";
  my $min = 1e300;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] < $min) { $min = $$vec[$i];  $ind = $i} }

  return $ind;
}

#---------------------------------------------------------
# vec_max
#---------------------------------------------------------
sub vec_max ($)
{
  my $vec = shift;

  my $max = -1e300;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] > $max) { $max = $$vec[$i]; } }

  return $max;
}

#---------------------------------------------------------
# vec_median
#---------------------------------------------------------
sub vec_median ($)
{
  my $vec = shift;

  my @full_vec;
  foreach my $num (@$vec)
  {
    if (length($num) > 0)
    {
      push(@full_vec, $num);
    }
  }

  @full_vec = sort { $a <=> $b } @full_vec;

  my $counts = @full_vec;

  # BEFORE WAS: return ($counts % 2 == 1) ? $full_vec[int($counts / 2)] : ($full_vec[int($counts / 2) - 1] + $full_vec[int($counts / 2)]) / 2
  # CHANGED NOV18,2007 because one property of median (e.g., comparing to mean) is that it is an original value of the set, and the above solution (for the even set size case) doesnt respect it...
  return $full_vec[int($counts / 2)];
}

#---------------------------------------------------------
# vec_true_median - here if the vec has even number of members
# we compute the mean of the two center members
#---------------------------------------------------------
sub vec_true_median ($)
{
  my $vec = shift;

  my @full_vec;
  foreach my $num (@$vec)
  {
    if (length($num) > 0)
    {
      push(@full_vec, $num);
    }
  }

  @full_vec = sort { $a <=> $b } @full_vec;

  my $counts = @full_vec;

  return ($counts % 2 == 1) ? ($full_vec[int($counts / 2)]) 
      : (($full_vec[int($counts / 2) - 1] + $full_vec[int($counts / 2)]) / 2);
}

#---------------------------------------------------------
# vec_quantile
#---------------------------------------------------------
sub vec_quantile ($$)
{
  my $vec = @_[0];
  my $quantile = @_[1];

  my @full_vec; 
  foreach my $num (@$vec)
  {
    if (length($num) > 0)
    {
      push(@full_vec, $num);
    }
  }

  @full_vec = sort { $a <=> $b } @full_vec;

  my $counts = @full_vec;

  my $save_quantile = $quantile;
  if ($quantile < 0)
  {
     $quantile = 0;
  }
  elsif ($quantile >= 1)
  {
     $save_quantile = (($counts - 1)/$counts);
  }

  return $full_vec[int($save_quantile * $counts)];
}

#---------------------------------------------------------
# vec_num_above_or_equal_min
#---------------------------------------------------------
sub vec_num_above_or_equal_min ($$)
{
    my $vec = @_[0];
    my $min = @_[1];

    my $res = 0;

    for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] >= $min) { $res++; } }

    return $res;
}

#---------------------------------------------------------
# vec_num_below_or_equal_max
#---------------------------------------------------------
sub vec_num_below_or_equal_max ($$)
{
    my $vec = @_[0];
    my $max = @_[1];

    my $res = 0;

    for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0 && $$vec[$i] <= $max) { $res++; } }

    return $res;
}

#---------------------------------------------------------
# mutual_information
#---------------------------------------------------------
sub vec_mutual_information ($$)
{
  my ($vecX, $vecY) = @_;



  DEBUG_VO("mutual_information between\n");
  DEBUG_VO("   X: @$vecX\n");
  DEBUG_VO("   Y: @$vecY\n");

  my $total = @$vecX;

  my %COUNTS_X;
  my %COUNTS_Y;
  my %COUNTS_XY;
  for (my $i = 0; $i < $total; $i++)
  {
         my $key_X = $$vecX[$i];
         my $key_Y = $$vecY[$i];
         my $key_XY = $$vecX[$i] . " " . $$vecY[$i];

         if (length($COUNTS_X{$key_X} == 0)) { $COUNTS_X{$key_X} = 1; } else { $COUNTS_X{$key_X} = $COUNTS_X{$key_X} + 1; }
         if (length($COUNTS_Y{$key_Y} == 0)) { $COUNTS_Y{$key_Y} = 1; } else { $COUNTS_Y{$key_Y} = $COUNTS_Y{$key_Y} + 1; }
         if (length($COUNTS_XY{$key_XY} == 0)) { $COUNTS_XY{$key_XY} = 1; } else { $COUNTS_XY{$key_XY} = $COUNTS_XY{$key_XY} + 1; }

         #DEBUG_VO("      $i: X[$key_X]=$COUNTS_X{$key_X}\n");
         #DEBUG_VO("         Y[$key_Y]=$COUNTS_Y{$key_Y}\n");
         #DEBUG_VO("         XY[$key_XY]=$COUNTS_XY{$key_XY}\n");
  }

  my $mi = 0;
  foreach my $key_XY (keys %COUNTS_XY)
  {
         my ($key_X, $key_Y) = split(" ", $key_XY);

         my $p_X = $COUNTS_X{$key_X} / $total;
         my $p_Y = $COUNTS_Y{$key_Y} / $total;
         my $p_XY = $COUNTS_XY{$key_XY} / $total;

         DEBUG_VO("      PX($key_X) = $p_X\n");
         DEBUG_VO("      PY($key_Y) = $p_Y\n");
         DEBUG_VO("      PXY($key_XY) = $p_XY\n");

         $mi += ($p_XY) * log($p_XY / ($p_X * $p_Y));
  }

  return $mi;
}

#---------------------------------------------------------
# vec_nth_statistic
#---------------------------------------------------------
sub vec_nth_statistic ($$)
{
  my $vec = @_[0];
  my $nth_statistic = @_[1];

  my @sorted_vec = sort { $a <=> $b } @$vec;

  return $sorted_vec[$nth_statistic];
}

#---------------------------------------------------------
# pearson
#---------------------------------------------------------
sub vec_pearson ($$)
{
  my ($vecX, $vecY) = @_;

  my $dot = 0;
  my $num = 0;
  for (my $i = 0; $i < @$vecX; $i++)
  {
    if (length($$vecX[$i]) > 0 and length($$vecY[$i]) > 0)
    {
      $dot += ($$vecX[$i] * $$vecY[$i]);
      $num++;
    }
  }

  my $pearson = $num > 0 ? $dot / $num : -1000;

  return $pearson;
}

#---------------------------------------------------------
# vec_stats
#---------------------------------------------------------
sub vec_permute (\@$)
{
  my ($vec_str) = @_;

  my @vec = @{$vec_str};

  my $n = @vec;

  while($n>1){
    $n--;
    my $new_index = int(rand($n+1));
    
    my $tmp = $vec[$n];
    $vec[$n] = $vec[$new_index];
    $vec[$new_index] = $tmp;
  }

  return @vec;
}

#---------------------------------------------------------
# vec_stats
#---------------------------------------------------------
sub vec_stats ($)
{
  my $vec = @_[0];

  my %stats = 0;
  for (my $i = 0; $i < @$vec; $i++) { if (length($stats{$$vec[$i]}) > 0) { $stats{$$vec[$i]} = $stats{$$vec[$i]} + 1; } else { $stats{$$vec[$i]} = 1; } }

  return %stats;
}

#---------------------------------------------------------
# vec_mad
#---------------------------------------------------------
sub vec_mad ($)
{
  my $vec = shift;

  my $true_med = vec_true_median($vec);
  # print "\n";
  my @diff_vec = ();
  for (my $i = 0; $i < @$vec; $i++) { push(@diff_vec, abs($$vec[$i] - $true_med));}
  #print "\n";
  #for (my $i = 0; $i < @$vec; $i++) { print "diff\t".$diff_vec[$i]."\n";}
  my $res = vec_true_median(\@diff_vec);

  return $res;
}

#---------------------------------------------------------
# vec_std
#---------------------------------------------------------
sub vec_std($)
{
  my $vec = shift;

  my $sum_x = 0;
  my $sum_xx = 0;
  for (my $i = 0; $i < @$vec; $i++) { $sum_x += $$vec[$i]; $sum_xx += $$vec[$i] * $$vec[$i]; }

  my $entries = vec_count_full_entries($vec);

  my $std = $entries > 0 ? $sum_xx / $entries - (($sum_x / $entries) * ($sum_x / $entries)) : "";

  return $std > 0 ? sqrt($std) : 0;
}


#---------------------------------------------------------
# vec_std
#---------------------------------------------------------
sub vec_stderr ($)
{
  my $vec = @_[0];

  my $sum_x = 0;
  my $sum_xx = 0;
  for (my $i = 0; $i < @$vec; $i++) { $sum_x += $$vec[$i]; $sum_xx += $$vec[$i] * $$vec[$i]; }

  my $entries = vec_count_full_entries($vec);

  my $stderr = $entries > 0 ? $sum_xx / $entries - (($sum_x / $entries) * ($sum_x / $entries)) : "";
  $stderr = $entries > 0 ? $stderr / $entries : $stderr;

  return sqrt($stderr);
}

#---------------------------------------------------------
# vec_sum
#---------------------------------------------------------
sub vec_sum($)
{
  my $vec = shift;

  my $sum = 0;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0) { $sum += $$vec[$i]; } }

  return $sum;
}

#---------------------------------------------------------
# vec_product
#---------------------------------------------------------
sub vec_product($)
{
  my $vec = shift;

  my $product = 1;
  for (my $i = 0; $i < @$vec; $i++) { if (length($$vec[$i]) > 0) { $product *= $$vec[$i]; } }

  return $product;
}

#---------------------------------------------------------
# vec_sum_log
#---------------------------------------------------------
sub vec_sum_log($)
{
  my $vec = shift;

  my $sum = 0;
  my $not_done = 1;
  my $first_val = -1;

  for (my $i = 0; ($i < @$vec) && ($not_done == 1); $i++)
  {
     if (length($$vec[$i]) > 0)
     {
	$sum = $$vec[$i];
	$first_val = $i;
	$not_done = 0;
     }
  }

  for (my $i = ($first_val + 1); $i < @$vec; $i++)
  {
     if (length($$vec[$i]) > 0)
     {
	my $x = $$vec[$i];
	my $y = $sum;
	if ($x >= $y)
	{
	   $y -= $x;
	}
	else
	{
	   my $t = $x;
	   $x = $y;
	   $y = $t - $x;
	}
	$sum = $x + log(1 + exp($y));
     }
  }

  return $sum;
}

#---------------------------------------------------------
# vec_union
#---------------------------------------------------------
sub vec_union ($$)
{
  my ($vecX, $vecY) = @_;

  my @res;
  my $counter = 0;

  my %h1;
  for (my $i = 0; $i < @$vecX; $i++) { $h1{$$vecX[$i]} = "1"; }
  for (my $i = 0; $i < @$vecY; $i++) { $h1{$$vecY[$i]} = "1"; }

  for my $key (keys %h1)
  {
    $res[$counter++] = $key;
  }

  return @res;
}

#---------------------------------------------------------
# vec_window_average
#---------------------------------------------------------
sub vec_window_average ($$)
{
  my $vec = @_[0];
  my $window_size = @_[1];

  my @res;

  for (my $i = 0; $i < @$vec; $i++)
  {
      my $start = $i - $window_size >= 0 ? $i - $window_size : 0;
      my $end = $i + $window_size < @$vec ? $i + $window_size : @$vec - 1;
      my $sum = 0;
      my $count = 0;
      for (my $j = $start; $j <= $end; $j++)
      {
	  if (length($$vec[$j])>0){
	    $sum += $$vec[$j];
	    $count++;
	  }
      }
      if ($count>0){
	$res[$i] = $sum / $count;
      }
      else{
	$res[$i] = "";
      }
  }

  return @res;
}

#---------------------------------------------------------
# vec_window_max
#---------------------------------------------------------
sub vec_window_max ($$)
{
  my $vec = @_[0];
  my $window_size = @_[1];

  my @res;

  for (my $i = 0; $i < @$vec; $i++)
  {
      my $start = $i - $window_size >= 0 ? $i - $window_size : 0;
      my $end = $i + $window_size < @$vec ? $i + $window_size : @$vec - 1;
      my $max = "";
      for (my $j = $start; $j <= $end; $j++)
      {
	  if (length($max) == 0 or $$vec[$j] > $max)
	  {
	      $max = $$vec[$j];
	  }
      }
      $res[$i] = $max;
  }

  return @res;
}

#---------------------------------------------------------
# vec_window_min
#---------------------------------------------------------
sub vec_window_min ($$)
{
  my $vec = @_[0];
  my $window_size = @_[1];

  my @res;

  for (my $i = 0; $i < @$vec; $i++)
  {
      my $start = $i - $window_size >= 0 ? $i - $window_size : 0;
      my $end = $i + $window_size < @$vec ? $i + $window_size : @$vec - 1;
      my $min = "";
      for (my $j = $start; $j <= $end; $j++)
      {
	  if (length($min) == 0 or $$vec[$j] < $min)
	  {
	      $min = $$vec[$j];
	  }
      }
      $res[$i] = $min;
  }

  return @res;
}

1
