#!/usr/bin/perl

use strict;
use Scalar::Util qw(looks_like_number);
use POSIX qw(floor ceil);

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/format_number.pl";
require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";


if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}

my $LOG10 = log(10);
my $LOG2 = log(2);

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

my $columns = get_arg("c", 0, \%args);
my $rows = get_arg("r", "All", \%args);
my $apply_to_all_columns = get_arg("A", 0, \%args);
my $min_filter = get_arg("min", "None", \%args);
my $max_filter = get_arg("max", "None", \%args);
my $mins_filter = get_arg("mins", "None", \%args);
my $maxs_filter = get_arg("maxs", "None", \%args);
my $abs_filter = get_arg("abs", "None", \%args);
my $below_abs_filter = get_arg("babs", "None", \%args);
my $str_filter = get_arg("str", "None", \%args);
my $equal_str_filter = get_arg("estr", "None", \%args);
my $equal_str_list_filter_str = get_arg("estr_list", "None", \%args);
my $not_equal_str_filter = get_arg("nstr", "None", \%args);
my $non_empty_filter = get_arg("ne", "None", \%args);
my $empty_filter = get_arg("empty", "None", \%args);
my $pairs_file = get_arg("pairs", "", \%args);
my $range_filter = get_arg("rng", "None", \%args);
my $range_strict_filter = get_arg("rngs", "None", \%args);
my $reverse_c = get_arg("reverse", 0, \%args);
my $inverse = get_arg("inv", 0, \%args);
my $add_number = get_arg("a", 0, \%args);
my $subtract_number = get_arg("s", 0, \%args);
my $multiply_number = get_arg("m", 1, \%args);
my $divide_number = get_arg("d", 1, \%args);

my $modulo_number = get_arg("mod", "", \%args);

my $absolute_value = get_arg("ab", "", \%args);
my $floor = get_arg("floor", 0, \%args);
my $ceil = get_arg("ceil", 0, \%args);
my $round = get_arg("round", 0, \%args);
my $sign = get_arg("sign", 0, \%args);
my $take_length = get_arg("len", 0, \%args);
my $take_log = get_arg("log", 0, \%args);
my $take_log2 = get_arg("log2", 0, \%args);
my $take_log10 = get_arg("log10", 0, \%args);
my $nonpos_log = get_arg("nonpl", "0", \%args);
my $take_exp_by_e = get_arg("e", 0, \%args);
my $take_exp = get_arg("exp", "", \%args);
my $take_pow = get_arg("pow", "", \%args);
my $divide_by_column = get_arg("dc", "", \%args);
my $multiply_by_column = get_arg("mc", "", \%args);
my $subtract_column = get_arg("sc", "", \%args);
my $negative_to_positive = get_arg("neg2pos", "", \%args);
my $add_column = get_arg("ac", "", \%args);
my $accumulate = get_arg("cumul", 0, \%args);
my $zero_div_value = get_arg("zero", "999999", \%args);
my $min_column = get_arg("minc", "", \%args);
my $max_column = get_arg("maxc", "", \%args);

my $set_min = get_arg("set_min", "", \%args);
my $set_max = get_arg("set_max", "", \%args);
my $skip = get_arg("skip", 0, \%args);
my $skip_columns = get_arg("skipc", 0, \%args);
my $substr = get_arg("substr", 0, \%args);
my $type = get_arg("type", "Number", \%args);
my $rm_regexp = get_arg("rmre", "", \%args);
my $to_regexp = get_arg("tore", "", \%args);
my $split_char = get_arg("splt_char", "", \%args);
my $pivot = get_arg("pv", 0, \%args);
my $pivot_val = get_arg("pv_val", 0.05, \%args);
my $pivot_high = get_arg("pv_high", 0, \%args);
my $pivot_low = get_arg("pv_low", 1, \%args);
my $split_del = get_arg("splt_d", "", \%args);
my $resolution = get_arg("res", "", \%args);
my $bins= get_arg("bins", "", \%args);
my $rescale= get_arg("rescale", "", \%args);
my $rescale_limits_set= get_arg("rescale_limits", "", \%args);
my $probdist= get_arg("pd", "", \%args);
my $subtract_mean= get_arg("z", "", \%args);
my $divide_mean= get_arg("dz", "", \%args);

my @accumulator;
my %columns_hash;
my @columns_list;

my $range_max;
my $range_min;
if ($range_filter ne "None"){
  ($range_min,$range_max)=($range_filter=~/^(\S+)\s*,\s*(\S+)/);
}
if ($range_strict_filter ne "None"){
  ($range_min,$range_max)=($range_strict_filter=~/^(\S+)\s*,\s*(\S+)/);
}

my %rows_hash;
my @rows_list = split(/\,/, $rows,-1);
foreach my $r (@rows_list) { $rows_hash{$r} = "1"; }

my $precision = get_arg("p", "", \%args);
my $add_str = get_arg("astr", "", \%args);
my $add_prefix_str = get_arg("bstr", "", \%args);
my $set_str = get_arg("set", "None", \%args);
my $quote = get_arg("quote", 0, \%args);

if ($quote)
{
	$add_prefix_str = "\"";
	$add_str = "\"";
}

my %pairs;
if (length($pairs_file) > 0)
{
  open(PAIRS_FILE, "<$pairs_file");
  while(<PAIRS_FILE>)
  {
    chop;

    my @row = split(/\t/,-1);

    $pairs{$row[0]}{$row[1]} = "1";
  }
}

my %equal_str_list_filter;
if (length($equal_str_list_filter_str) > 0)
{
  my @row = split(/;/, $equal_str_list_filter_str);
  for (my $i = 0; $i < @row; $i++)
  {
    $equal_str_list_filter{$row[$i]} = "1";
  }
}

my @headers;


my @col_max_array;
my @col_min_array;
my @col_sum_array;
my @col_possum_array;
my @col_count_array;
my $rescale_top;
my $rescale_bottom;
my $tmp_id=int(rand(100000000));
my $prescan_file=0;
if ($bins ne "" or $rescale ne "" or $probdist ne "" or $subtract_mean ne "" or $divide_mean ne ""){ $prescan_file=1 }
my $rescale_limits_set_max ;
my $rescale_limits_set_min ;
if ($rescale ne ""){
  ($rescale_bottom, $rescale_top)=(split /\s*,\s*/,$rescale,-1);
}
if ($rescale_limits_set ne ""){
  $prescan_file=0;
  ($rescale_limits_set_min, $rescale_limits_set_max)=(split /\s*,\s*/,$rescale_limits_set,-1);
}
if ($skip >= 1)
{
  my $line = <$file_ref>;
  print "$line";
  chomp $line;
  @headers = split(/\t/, $line,-1);
}

for (my $i = 1; $i < $skip; $i++)
{
  my $line = <$file_ref>;
  print "$line";
}


if ($prescan_file){
  open (TMP,">tmp_$tmp_id");
  while (my $r=<$file_ref>){
    print TMP $r;
    chomp $r;
    my @row=split /\t/,$r,-1;
    if ($#columns_list==-1){
      @columns_list=parseRanges($columns,$#row);
      foreach my $c (@columns_list) { $columns_hash{$c} = "1"; }
    }
    for (my $i = 0; $i < scalar(@row); $i++) {
      if (($columns_hash{$i} eq "1" or $apply_to_all_columns == 1) and pass_filter($row[0], $headers[$i], $row[$i])){
	$col_sum_array[$i]+=$row[$i];
	$col_count_array[$i]++;
	if ($row[$i]>0){
	  $col_possum_array[$i]+=$row[$i];
	}
	if (!defined($col_max_array[$i]) or $col_max_array[$i]<$row[$i]){
	  $col_max_array[$i]=$row[$i];
	}
	if (!defined($col_min_array[$i]) or $col_min_array[$i]>$row[$i]){
	  $col_min_array[$i]=$row[$i];
	}
      }
    }
  }
  close (TMP);
  open (TMP,"tmp_$tmp_id");
  $file_ref=\*TMP;
}

my $row_count=$skip-1;

while (my $line=<$file_ref>)
{
  chomp $line;
  $row_count++;

  if (!$rows_hash{All} and !$rows_hash{$row_count}) {
    print "$line\n";
    next;
  }

  my @row = split /\t/,$line,-1;
  if ($#columns_list==-1){
    @columns_list=parseRanges($columns,$#row);
    foreach my $c (@columns_list) { $columns_hash{$c} = "1"; }
  }


  for (my $i = 0; $i < $skip_columns; $i++)
  {
    print "$row[$i]\t";
  }


  for (my $i = $skip_columns; $i < scalar(@row); $i++)
  {
    if (($columns_hash{$i} eq "1" or $apply_to_all_columns == 1) and &pass_filter($row[0], $headers[$i], $row[$i]))
    {
      my $num = $row[$i];

      if ($substr != 0 )
      {
	  if (length($num) > $substr)
	  {
	      $num = substr($num, 0, $substr);
	  }
      }
      elsif (length($divide_by_column) > 0)
      {
	if (length($num) > 0)
	{ 
	  if ($row[$divide_by_column] != 0) { $num /= $row[$divide_by_column]; }
	  elsif ($num == 0) { $num = 0; }
	  else { $num = $zero_div_value; }
	}
      }
      elsif (length($multiply_by_column) > 0)
      {
	if (length($num) > 0) { $num *= $row[$multiply_by_column]; }
      }
      elsif ($inverse)
      {
	if (length($num) > 0) {
	  if ($num==0){ $num=$zero_div_value }
	  else{ $num=1/$num }
	}
      }
      elsif (length($rescale) > 0)
      {
	if (length($num) > 0) {
	  if ($rescale_limits_set ne ""){
	    $col_max_array[$i]=$rescale_limits_set_max;
	    $col_min_array[$i]=$rescale_limits_set_min;
	  }
	  if ($col_max_array[$i]==$col_min_array[$i]){
	    $num=$rescale_bottom;
	  }
	  else{
	      $num = ( $num-$col_min_array[$i]) * ( ($rescale_top-$rescale_bottom)/($col_max_array[$i]-$col_min_array[$i]) ) + $rescale_bottom;
	    }
	}
      }
      elsif (length($bins) > 0)
      {
	if (length($num) > 0) {
	  my $binsize=($col_max_array[$i]-$col_min_array[$i])/$bins;
	  my $f = 0.5;
	  if ($num==$col_max_array[$i]){
	    $f = -0.5
	  }
	  $num = (int( ($num-$col_min_array[$i])/$binsize )+$f) * $binsize + $col_min_array[$i];
	}
      }
      elsif (length($probdist) > 0)
      {
	if (length($num) > 0) {
	  if ($num<0) {
	    $num=0;
	  }
	  else{
	    $num = ($col_possum_array[$i] == 0 ? 0 : $num/$col_possum_array[$i]);
	  }
	}
      }
      elsif (length($subtract_mean) > 0)
      {
		if (length($num) > 0) {
			$num -=  $col_sum_array[$i]/$col_count_array[$i];
		}
      }
      elsif (length($divide_mean) > 0)
      {
		if (length($num) > 0) {
			$num /=  $col_sum_array[$i]/$col_count_array[$i];
		}
      }
      elsif (length($subtract_column) > 0)
      {
	  if (length($num) > 0) { $num -= $row[$subtract_column]; }
      }
      elsif (length($add_column) > 0)
      {
	  if (length($num) > 0) { $num += $row[$add_column]; }
      }
      elsif (length($min_column) > 0)
      {
	  $num = ($num < $row[$min_column]) ? $num : $row[$min_column];
      }
      elsif (length($max_column) > 0)
      {
	  $num = ($num < $row[$max_column]) ? $row[$max_column] : $num; 
      }
      elsif ($take_length){
	$num = length($num);
      }
      elsif ($reverse_c) {
      $num = reverse ($num);
      }
      elsif ($pivot == 1)
      {
	if (length($num) > 0 && $num =~ /^[0-9\.-]*[Ee]?[+0-9\.-]*$/)
	{
	  if ($num > $pivot_val)
	  {
	    $num = $pivot_high;
	  }
	  else
	  {
	    $num = $pivot_low;
	  }
	}
      }
      elsif (length($num) > 0 and $type eq "Number" and length($add_str) == 0 and length($add_prefix_str) == 0 and $set_str eq "None" && length($rm_regexp) == 0 && length($split_del) == 0 && length($split_char) == 0)
      {
	if ($accumulate){ $accumulator[$i] += $num; $num = $accumulator[$i] }
	$num += $add_number;
	$num -= $subtract_number;
	$num *= $multiply_number;
	$num = ($divide_number == 0 ? $zero_div_value : $num/$divide_number);
	if ($modulo_number) { $num %= $modulo_number ; }
	if ($take_log == 1) { $num = $num > 0 ? log($num) : $nonpos_log; }
	if ($take_log2 == 1) { $num = $num > 0 ? log($num) / $LOG2 : $nonpos_log; }
	if ($take_log10 == 1) { $num = $num > 0 ? log($num) / $LOG10 : $nonpos_log; }
	if ($absolute_value ne "" and $num<0) { $num=$num*(-1); }
	if ($floor == 1) { $num = floor($num); }
	if ($ceil == 1) { $num = ceil($num); }
#	if ($round == 1) { $num = int($num) + (($num - int($num)))/0.5 >=1 ? 1 : 0); }
	if ($round == 1) { $num = abs(floor($num) - $num) < abs(ceil($num) - $num) ? floor($num) : ceil($num); }
	if ($sign == 1) { $num = $num == 0 ? 0 : $num / abs($num); }
	if ($take_exp_by_e == 1) { $num = exp($num); }
	if (length($take_exp) > 0) { $num = $take_exp ** $num; }
	if (length($take_pow) > 0) { $num = $num == 0 ? 0 : $num ** $take_pow; }
	if (length($set_min) > 0 and $num < $set_min) { $num = $set_min; }
	if (length($set_max) > 0 and $num > $set_max) { $num = $set_max; }
	if (length($negative_to_positive) > 0 and $num < 0) { $num = -$num; }
	if (length($resolution) > 0) {
	  my $tmp = $resolution * &format_number($num / $resolution, 0); $num = $num - $tmp < $tmp + $resolution - $num ? $tmp : $tmp + $resolution; }
      }
      elsif (length($rm_regexp) > 0)
      {
        $num =~ s/$rm_regexp/$to_regexp/g;
      }
      elsif (length($split_char) > 0)
      {
		my @split_col = split (//, $num);
	 	$num = join ("$split_char", @split_col);
      }
      elsif (length($split_del) > 0)
      {
	 my @split_col = split (/$split_del/, $num,-1);
	 $num = join ("\t", @split_col);
      }
      else
      {
	if ($set_str ne "None")
	{
	  if ($set_str eq "EMPTY") { $num = ""; }
	  else { $num = $set_str; }
	}
	else
	{
	  $num = "$add_prefix_str$num$add_str";
	}
      }

      if ( length($precision) > 0 and looks_like_number($num) ) { $num = &format_number($num, $precision); }

      print "$num";

      #print STDERR "MODIFYING [$row[0], $headers[$i]] --> VALUE [$row[$i]] --> [$num]\n";
    }
    else
    {
      print "$row[$i]";
    }

    if ($i < scalar(@row) - 1) { print "\t"; }
  }

  print "\n";
}

sub pass_filter
{
  my ($row_name, $column_name, $num) = @_;

  my $pass = 1;

  my $sci_number = $num =~ /^[\-]?[0-9\.]+[Ee][\-]?[0-9]+/;

  if ($num =~ /[A-Z]/ and $sci_number != 1) { if ($min_filter ne "None" or $max_filter ne "None" or $mins_filter ne "None" or $maxs_filter ne "None" or $range_filter ne "None" or $range_strict_filter ne "None") { $pass = 0; } }

  if ($min_filter ne "None" and $num < $min_filter) { $pass = 0; }
  if ($max_filter ne "None" and $num > $max_filter) { $pass = 0; }
  if ($mins_filter ne "None" and $num <= $mins_filter) { $pass = 0; }
  if ($maxs_filter ne "None" and $num >= $maxs_filter) { $pass = 0; }
  if ($range_filter ne "None" and ($range_max<$num or $range_min>$num)) { $pass = 0; }
  if ($range_strict_filter ne "None" and ($range_max<=$num or $range_min>=$num)) { $pass = 0; }
  if ($abs_filter ne "None" and abs($num) < $abs_filter) { $pass = 0; }
  if ($below_abs_filter ne "None" and abs($num) > $below_abs_filter) { $pass = 0; }
  if ($str_filter ne "None" and not($num =~ /$str_filter/)) { $pass = 0; }
  if ($equal_str_filter ne "None" and $num ne $equal_str_filter) { $pass = 0; }
  if ($equal_str_list_filter_str ne "None" and $equal_str_list_filter{$num} ne "1") { $pass = 0; }
  if ($not_equal_str_filter ne "None" and $num eq $not_equal_str_filter) { $pass = 0; }
  if ($non_empty_filter ne "None" and length($num) == 0) { $pass = 0; }
  if ($empty_filter ne "None" and length($num) > 0) { $pass = 0; }
  if (length($pairs_file) > 0 and $pairs{$row_name}{$column_name} ne "1") { $pass = 0; }

  return $pass;
}

if ($prescan_file){
  close TMP;
  unlink "tmp_$tmp_id";
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/modify_column.pl <file>

   Modifies a column according to predefined operations

   WHAT TO MODIFY:

      -c <num>       The column to modify (default: 0; specify multiple columns using commas) (zero-based)
      -r <num>       The rows to modify (default: All except -skipped; specify multiple rows using commas) (zero-based)
      -A             Apply transformation to ALL columns
      -min <num>     Apply transformation only if entry is above (or equal) num
      -max <num>     Apply transformation only if entry is below (or equal) num
      -mins <num>    Apply transformation only if entry is strictly above num
      -maxs <num>    Apply transformation only if entry is strictly below num
      -rng <num1,num2> Apply transformation only if num1<=entry<=num2
      -rngs <num1,num2> Apply transformation only if num1<entry<num2
      -abs <num>:    Apply transformation only if entry is above <num> or below -<num>
      -babs <num>:   Apply transformation only if entry is above -<num> and below <num>
      -str <str>:    Apply transformation only if entry contains <str>
      -estr <str>:   Apply transformation only if entry is equal to <str>
      -nstr <str>:   Apply transformation only if entry is *not* equal to <str>
      -ne:           Apply transormation only if entry is not empty
      -empty:        Apply transormation only if entry is empty
      -pairs <file>: Apply transformation only if the entry is from a row and column 
                     that appear in a line of <file> as <row name> <tab> <column name> 
      -skip <num>:   Skip num rows in the file (default: 0)
      -skipc <num>:  Skip num columns (when applying to all columns, can skip the first few) (default: 0)

   HOW TO MODIFY:

      -a <num>       Add <num> to the column (default: 0)
      -s <num>       Substract <num> to the column (default: 0)
      -m <num>       Multiply column by <num> (default: 1) (write '"-1"' to pass negative numbers)
      -d <num>       Divide column by <num> (default: 1) (write '"-1"' to pass negative numbers)
      -mod <num>     Modulo column by <num> (default: 1) (write '"-1"' to pass negative numbers)
      -inv           Inverse of column
      -ab            Absolute value of column
      -floor         Floor of column values
      -ceil          Ceiling of column values
      -round         Round to the nearest integer. Numbers that are halfway between two integers are rounded down (like int())
      -sign          Sign of column values (+1 for positive, -1 for negative, 0 for zero)
      -neg2pos       Convert negative numbers to positive
      -log           Take log (natural base) of the column
      -log2          Take log (base 2) of the column
      -log10         Take log (base 10) of the column
      -nonpl <str>   String to print out in case of non-positive log (default: "0")
      -exp <num>     Take num ** (value of the column)
      -e             Take e=2.1782... ** (value of column)
      -pow <num>     Raise vaule to the power of num (value of the column)
      -substr <num>  If the column has more than num chars, then substr the first num
      -set_min <num> If the column has a value less than num, set it to num
      -set_max <num> If the column has a value greater than num, set it to num
      -type <str>:   <Number/String> (default: Number)
      -len           Replace value (treated as a string) by its length in characters
      -reverse       Reverse the string (ABC becomes CBA)

      -dc <col>      Divide the column by the number in column <col>
      -mc <col>      Multiply the column by the number in column <col>
      -ac <col>      Add the number in column <col> to the column
      -sc <col>      Subtract the number in column <col> from the column
      -minc <col>    Take the min of the column and <col>
      -maxc <col>    Take the max of the column and <col>

      -bins <num>:   discretize the column into <num> bins.
      -rescale <str> given as 'bottom,top'. linearly rescales the data in the column so that the smallest value will
                     be <bottom> and the largest will be <top>. If the data is equal in all entries, all entries
                     will be set to <bottom>.
      -rescale_limits <str>: given as 'min,max'. when used with -rescale assumes that the maximal value of the column is
                     max and the minimal value is min.
      -pd            makes the column a probability distribution by changing negative values to zero and
                     then dividing by the sum of the column.
      -z             Subtract mean
      -dz            Divide by mean

      -p <num>       Precision (take only <num> sig. digits. default: don't fix)
      -res <num>     Convert numbers to a <num> resolution (e.g., 0.05)
      -astr <str>    Add <str> to the end of the column (default: "")
      -bstr <str>    Add <str> to the beginning of the column (default: "")
      -quote         Quote the column (add quotes before and after string)
      -set <str>     Set the entry to <str> (put EMPTY for getting an empty entry)
      -rmre <regexp> Remove string that matches the regular expression <regexp> from column

      -pv <bool>     Set numerical values around a given pivot (see -pv_val, -pv_high, -pv_low). Activated when <bool> = 1.
                     Non numerical values are not modified.
      -pv_val  <num> The pivot value (default: 0.05)     
      -pv_high <num> The number to set values greater then the pivot value (default: 0)
      -pv_low  <num> The number to set values smaller equal then the pivot value (default: 1)

      -splt_d <del>  Split the column by the delimiter <del> (default: none).
      -cumul         Accumulate (row i = row 1 + .. + row i)
      -zero <num>    When dividing by zero, uses <num> instead of giving an error (default: 999999)


