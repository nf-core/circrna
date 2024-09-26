#!/usr/bin/perl

use strict;
use POSIX;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/ltqnorm.pl";

my $MAXDOUBLE = 1e200;

my $LOG10     = log(10);

my $PI = 3.141592654;

my $MIN_NUMBER_TO_EXPONENTIATE = -100;
   
#print ComputeHyperPValue($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]) . "\n";
#print NormalStd2Pvalue($ARGV[0]);

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeBinomial
{
   my ($p, $n, $r) = @_;

   my $binom = -$MAXDOUBLE;
   my $x=1;
   for (my $i = $r; $i <= $n; $i++)
   {
     $binom = &AddLog($binom, lchoose($i, $n) + $i * log($p) + ($n - $i) * log(1-$p));
   }

   return exp($binom);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeBinomial2
{
   my ($p, $n, $r) = @_;

   my $binom = -$MAXDOUBLE;
   my $x=1;
   for (my $i = $r; $i <= $n; $i++)
   {
     $binom = &AddLog($binom, lchoose_Stanica($i, $n) + $i * log($p) + ($n - $i) * log(1-$p));
   }

   return exp($binom);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Log
{
   my ($x) = @_;

   if ($x == 0.0)
   {
      return -$MAXDOUBLE;
   }
   else
   {
      return log($x);
   }
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeTTest
{
   my ($x, $dof) = @_;

   my $u = $dof / ($dof + $x * $x);

   my $val = 0.5 * &BetaCDF($u, 0.5 * $dof, 0.5);

   return $x > 0 ? 1 - $val : $val;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeTValue
{
   my ($n1, $m1, $s1, $n2, $m2, $s2) = @_;

   my $se = sqrt ((($s1 * $s1) / $n1) + (($s2 * $s2) / $n2));

   my $val = ($m1 - $m2) / $se;

   return $val < 0 ? $val : -$val;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub BetaCDF
{
   my ($x, $a, $b) = @_;
  
   my $bt;

   if ($x == 0 or $x == 1) { $bt = 0; }
   else
   {
     $bt = exp(&lgamma($a + $b) - &lgamma($a) - &lgamma($b) + $a * log($x) + $b * log(1 - $x)); 
   }

   return $x < (($a + 1) / ($a + $b + 2)) ? return ($bt * &BetaCF($x, $a, $b) / $a) : 1 - ($bt * &BetaCF(1 - $x, $a, $b) / $b)
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub BetaCF
{
   my ($x, $a, $b) = @_;

   my $maxit = 100;

   my $eps = 0.0000003;
   my $am = 1;
   my $bm = 1;
   my $az = 1;
   my $qab = $a + $b;
   my $qap = $a + 1;
   my $qam = $a - 1;
   my $bz = 1 - $qab * $x / $qap;
   my $tem;
   my $em;
   my $d;
   my $bpp;
   my $bp;
   my $app;
   my $aOld;
   my $ap;

   for (my $m = 1; $m <= $maxit; $m++)
   {
     $em = $m;
     $tem = $em + $em;
     $d = $em * ($b - $m) * $x / (($qam + $tem) * ($a + $tem));
     $ap = $az + $d * $am;
     $bp = $bz + $d * $bm;
     $d = -($a + $em) *($qab + $em) * $x / (($a + $tem) * ($qap + $tem));
     $app = $ap + $d * $az;
     $bpp = $bp + $d * $bz;
     $aOld = $az;
     $am = $ap / $bpp;
     $bm = $bp / $bpp;
     $az = $app / $bpp;
     $bz = 1;
     if (&Abs($az - $aOld) < $eps * &Abs($az))
     {
       last; 
     }
   }

   return $az;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Abs
{
  my ($x) = @_;

  return $x > 0 ? $x : -$x; 
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeLog10HyperPValue
{
   my ($k, $n, $K, $N) = @_;

   my $p = $K / $N;

   if ($k <  $p * $n)
   {
      return 0;
   }
   
   
   if ($N == 0 or $n == 0 or ($k > $n) or ($K > $N))
   {
      return undef;
   }

   my $pvalue = &lchoose($k, $K) + &lchoose($n - $k, $N - $K) - &lchoose($n, $N);

   return $pvalue / $LOG10;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeHyperPValue
{
   my ($k, $n, $K, $N) = @_;

   my $p = $K / $N;

   my $d = ($k <  $p * $n) ? -1 : +1;

   my $PVal = 0;

   if ($d == 1)
   {
       $PVal = -$MAXDOUBLE;

       for (; $k >= 0 and $k <= $n and $k <= $K ; $k += $d)
       {
	   my $x = &lchoose($k, $K) + &lchoose($n - $k, $N - $K) - &lchoose($n, $N);

#	   my $x = &log_bincoef_approx($k, $K) + &log_bincoef_approx($n - $k, $N - $K) - &log_bincoef_approx($n, $N);

	   $PVal = &AddLog($PVal, $x);
       }
   }

   return exp($PVal);
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeHyperPValue2
{
  my ($k, $n, $K, $N) = @_;
  my $p = $K / $N;
  my $d = ($k <  $p * $n) ? -1 : +1;
  my $PVal = 0;
   if ($d == 1)
     {
       $PVal = -$MAXDOUBLE;
       for (; $k >= 0 and $k<=$n and $k<=$K; $k += $d)
       {
           my $x = &lchoose_Stanica($k, $K) + &lchoose_Stanica($n - $k, $N - $K) - &lchoose_Stanica($n, $N);
           $PVal = &AddLog($PVal, $x);
         }
     }
   return exp($PVal);
}

# bound formulas taken from "Good lower and upper bounds on binomial coefficients" JIPAM (2001) by Stanica
sub lchoose_Stanica{
  my $n=shift;
  my $k=shift;
  return (lchoose_Stanica_upperb($n,$k)+lchoose_Stanica_lowerb($n,$k))/2;
}

#upper bound
sub lchoose_Stanica_upperb{
  my $n=shift;
  my $k=shift;
  if ($k==$n or $n==0){return 0};
  my $m=$k/$n;
  return log(1/sqrt(2*(4*atan2(1,1))))-0.5*log($n)+($m*$n+0.5)*log($m)-(($m-1)*$n+0.5)*log($m-1)-($n+0.5)*log(1) ;
}

#lower bound
sub lchoose_Stanica_lowerb{
  my $n=shift;
  my $k=shift;
  if ($k==$n or $n==0){return 0};
  my $m=$k/$n;
  return log(1/sqrt(2*(4*atan2(1,1))))+1-1/(8*$n)-0.5*log($n)+($m*($n-1)+1)*log($m)-($m-1)*($n-1)*log($m-1)-($n+0.5)*log(1) ;
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub log_bincoef_approx
  {
    my ($k,$n);

    my $MY_PI = 4 * atan2(1, 1);

    if ($k==$n or $k==0){
      return 0;
    }

    my $m = $n/$k;

    my $const = (1 / sqrt(2 * $MY_PI));

    return log($const) + log($k)*(-0.5) + log($m)*($m*($k-1)+1) - log($m-1)*($m-1)*($k-1);

}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub log10 {
  my $n = shift;
  return log($n)/log(10);
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeHyperPValueOneTailed
{
   my ($k, $n, $K, $N) = @_;

   my $p = $K / $N;

   my $d = ($k <  $p * $n) ? -1 : +1;

   my $PVal = 0;

   if ($d == 1)
   {
       $PVal = -$MAXDOUBLE;

       for (; $k >= 0 and $k <= $n; $k += $d)
       {
	   my $x = &lchoose($k, $K) + &lchoose($n - $k, $N - $K) - &lchoose($n, $N);

	   $PVal = &AddLog($PVal, $x);
       }
   }

   return exp($PVal);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeWilcoxonSumRankTest
{
   my ($m, $n, $w) = @_;

   # m = number of items in sample 1
   # n = number of items in sample 2
   # w = sum of the ranks in sample 2

   my $mean = $n * ($m + $n + 1) / 2;
   my $std = sqrt($m * $n * ($m + $n + 1) / 12);
   
   my $num_std = abs(($w - $mean) / $std);
	
	#print "m=$m n=$n w=$w $mean\t$std\t$num_std\n";

   return &NormalStd2Pvalue($num_std);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeMannWhitneyUTest
{
	my ($n1, $n2, $r1) = @_;
	
	my $U = $n1*$n2 + ($n1*($n1+1)/2) - $r1;
	return ($U / ($n1*$n2));
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ComputeWilcoxonSumRankTestOneTailed
{
   my ($m, $n, $w) = @_;

   # m = number of items in sample 1
   # n = number of items in sample 2
   # w = sum of the ranks in sample 2
   #
   # Checking hypothesis that sample 1 has sum of ranks lower than
   # sample 2.

   my $mean = $n * ($m + $n + 1) / 2;
   my $std = sqrt($m * $n * ($m + $n + 1) / 12);
   
   my $num_std = abs(($w - $mean) / $std);
	
	#print "m=$m n=$n w=$w $mean\t$std\t$num_std\n";

   my $twoTailed = &NormalStd2Pvalue($num_std);
   
   return ($w > $mean ? $twoTailed : (1-$twoTailed));
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub lchoose
{
  my ($k, $N) = @_;

  return $k > $N ? 0 : &lgamma($N + 1) - &lgamma($k + 1) - &lgamma($N - $k + 1);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub AddLog
{
   my ($x, $y) = @_;

   if ($x == -$MAXDOUBLE) { return $y };
   if ($y == -$MAXDOUBLE) { return $x };
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

   return $x + log(1 + exp($y));
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub lgamma
{
  # per code from numerical recipies
  my $xx = $_[0];

  my ($j, $ser, $stp, $tmp, $x, $y);
  my @cof = (0.0, 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5);
  my $stp = 2.5066282746310005;

  $x = $xx;
  $y = $x;
  $tmp = $x + 5.5;
  $tmp = ($x+0.5)*log($tmp)-$tmp;
  $ser = 1.000000000190015;
  foreach $j ( 1 .. 6 )
  {
    $y+=1.0;
    $ser+=$cof[$j]/$y;
  }
  return $tmp + log($stp*$ser/$x);
}

#-----------------------------------------------------------------------------
# $double ComputeSymmetricUniformCdf(\@list reals, $int already_sorted=0)
#-----------------------------------------------------------------------------
sub ComputeSymmetricUniformCdf
{
   my ($reals, $already_sorted) = @_;
   $already_sorted = defined($already_sorted) ? $already_sorted : 0;

   my $list;
   if(not($already_sorted))
   {
      my @list = sort {$b <=> $a} @{$reals};
      push(@list, 0.0);
      $list = \@list;
   }
   else
   {
      $list = $reals;
   }

   if(scalar(@{$list}) == 1)
   {
      return 1.0;
   }
   else
   {
      my $cummulative_probability = 0.0;

      for(my $i = 0; $i < scalar(@{$list}) - 1; $i++)
      {
         my $range = $$list[$i] - $$list[$i + 1];

         my @list = @{$list};

         splice(@list, $i, 1);

         $cummulative_probability += $range * &ComputeSymmetricUniformCdf(\@list, 1);
      }

      return $cummulative_probability;
   }
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Pearson2FisherZ
{
   my ($r) = @_;

   my $z = defined($r) ? 0.5 * (log(1 + $r) - log(1 - $r)) : undef;

   return $z;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Pearson2FisherZscore
{
   my ($r, $dimensions) = @_;

   my $z      = &Pearson2FisherZ($r);
   my $zscore = (defined($z) and $dimensions > 3) ? $z * sqrt($dimensions - 3) : undef;

   return $zscore;
}

my $TWICE_A = 1.7155277699214135;

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub sample_normal
{
  while (1)
  {
    my $u = rand();
    if ($u != 0.0)
    {
      my $v = (rand() - 0.5) * $TWICE_A;
      my $x = $v / $u;
      my $sqr_x = $x*$x;
      if ($sqr_x <= 6 - 8*$u + 2*$u*$u)
      {
        return $x;
      }
      if (!($sqr_x >= 2 / $u - 2 * $u))
      {
        if ($sqr_x <= -4 * log($u))
        {
          return $x;
        }
      }
    }
  }
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub sample_normal_general
{
  my ($mean, $std) = @_;

  my $normal = &sample_normal();

  return ($std * $normal) + $mean;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub compute_normal_general
{
  my ($normal, $mean, $std) = @_;

  return (1 / (sqrt(2 * $PI) * $std)) * exp(-($normal - $mean) * ($normal - $mean) / (2 * $std * $std));
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub compute_log_normal_general
{
  my ($normal, $mean, $std) = @_;

  my $prob = &compute_normal_general($normal, $mean, $std);

  return $prob == 0 ? -$MAXDOUBLE : log($prob);
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub NormalStd2Pvalue
{
  my ($std) = @_;

  my $max = 1;
  my $min = 0;
  my $res = ($min + $max) / 2;

  while ($max - $min > 0.000000000000001)
  {
    my $qnorm = &ltqnorm($res);

    if ($qnorm < $std)
    {
      $min = $res;
    }
    else
    {
      $max = $res;
    }

    $res = ($max + $min) / 2;
  }

  return 1 - $res;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub IntersectNormals
{
   my ($mean1, $var1, $mean2, $var2, $lower_solution) = @_;

   my $a = $var1 - $var2;
   my $b = 2 * $mean1 * $var2 - 2 * $mean2 * $var1;
   my $c = $mean2 * $mean2 * $var1 - $mean1 * $mean1 * $var2 + $var1 * $var2 * 0.5 * log($var2 / $var1);

   my $res = 0;

   if ($a == 0)
   {
       $res = -$c / $b;
   }
   else
   {
       my $solution1 = (-$b + sqrt($b * $b - 4 * $a * $c)) / (2 * $a);
       my $solution2 = (-$b - sqrt($b * $b - 4 * $a * $c)) / (2 * $a);

       $res = $lower_solution == 1 ? &Min($solution1, $solution2) : &Max($solution1, $solution2);
   }

   return $res;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Min
{
    my ($num1, $num2) = @_;

    return $num1 < $num2 ? $num1 : $num2;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Max
{
    my ($num1, $num2) = @_;

    return $num1 > $num2 ? $num1 : $num2;
}

#-----------------------------------------------------------------------------
# compute log(exp(a_1)+exp(a_2)+...exp(a_n)) using:
# max(a_1,a_2,..a_n) + log(1+exp(a_2-max(a_1,a_2,..a_n))+...exp(a_n-max(a_1,a_2,..a_n)))
#-----------------------------------------------------------------------------
sub LogOfSumOfExps 
{
    my ($ref_exps) = @_;

    my @exps = @$ref_exps;
    my $num_exps = scalar(@exps);
    my $max_index = 0;
    my $max = $exps[0];
    my $i = 0;
    
    for ($i=1; $i<$num_exps; $i++)
    {
	if ($exps[$i] > $max)
	{
	    $max = $exps[$i];
	    $max_index = $i;
	}
    }

    my $remaining_summand = 1.0;

    for ($i = 0; $i < $num_exps; $i++)
    {
	if ($i != $max_index)
	{
	    my $difference = $exps[$i] - $max;
	    if ($difference > $MIN_NUMBER_TO_EXPONENTIATE){
		$remaining_summand += exp($difference);
	    }
	}
    }
    
    my $result = $max + log($remaining_summand);
    
    return $result;
}

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub KolmogorovSmirnovProbability
{
   my ($size1, $size2, $res_d_statistic) = @_;

   my $dtmp1 = sqrt(($size1 * $size2)/($size1 + $size2));

   my $x = ($dtmp1 + 0.12 + 0.11 / $dtmp1) * $res_d_statistic;

   my $x2 = -2.0 * $x * $x;
   my $eps1 = 10 * DBL_EPSILON;
   my $eps2 = DBL_EPSILON;

   my $i = 1;
   my $fac = 2.0;
   my $sum = 0.0;
   my $term = 0.0;
   my $absterm = 0.0;

   do
   {
      $absterm = abs($term);
      $term = $fac * exp($x2 * $i * $i);
      $sum += $term;
      $fac = -1 * $fac;
      $i++;
   }
   while (abs($term) > $eps1 * $absterm && abs($term) > $eps2 * $sum);

   return $sum;

}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ArrayRefMean {
  my $array_ref = shift;
  my $num_elements = scalar @$array_ref;
  my $sum = 0;

  foreach (@$array_ref) {
    $sum += $_;
  }
  return $sum / $num_elements;
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub ArrayRefStd {
  my $array_ref = shift;
  my $num_elements = scalar @$array_ref;
  my $sum = 0;
  my $sumsq = 0;

  foreach (@$array_ref) {
    $sum += $_;
    $sumsq += ($_ **2);
  }
  return sqrt( $sumsq/$num_elements -
               (($sum/$num_elements) ** 2));
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Bonferroni {
  my ($alpha,$num_pvalues) = @_;
  if ( $num_pvalues == 0 ) { return $alpha; }
  return $alpha/$num_pvalues;
}


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
sub Fdr {
  my ($alpha,$num_pvalues,$pvalues_ref) = @_;
  if ( $num_pvalues == 0 ) { return 0; }

  my @sorted_pvalues = sort { $a <=> $b } @$pvalues_ref;

  for ( my $i = $num_pvalues-1 ; $i >= 0 ; $i-- ) {
    my $randomly_expected = $sorted_pvalues[$i] * $num_pvalues;
    my $current_pvalue = $randomly_expected / ($i+1);
    if ( $current_pvalue <= $alpha ) { return $sorted_pvalues[$i]; }
  }
  return 0;
}


    
#print NormalStd2Pvalue($ARGV[0]) . "\n";
#print ComputeWilcoxonSumRankTest($ARGV[0], $ARGV[1], $ARGV[2]) . "\n";
#print ComputeBinomial($ARGV[0], $ARGV[1], $ARGV[2]) . "\n";
#print ComputeHyperPValue($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]) . "\n";  #(k, n, K, N)
#print ComputeLog10HyperPValue($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]) . "\n";  #(k, n, K, N)
#print &ComputeTTest($ARGV[0], $ARGV[1]) . "\n";
#print &compute_normal_general($ARGV[0], $ARGV[1], $ARGV[2]) . "\n";
#print &compute_log_normal_general($ARGV[0], $ARGV[1], $ARGV[2]) . "\n";
#print &IntersectNormals(0, 1, 2, 4, 1);
#print "\n";
#for(my$i=0;$i<=256;$i++){print "$i\t" . chr($i) . "A\n";}

1

