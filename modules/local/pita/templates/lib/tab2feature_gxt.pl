#!/usr/bin/perl

use strict;

require "/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/load_args.pl";

if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}

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

my $name = get_arg("n", "Gene names", \%args);
my $desc = get_arg("d", "", \%args);
my $org = get_arg("o", "", \%args);
my $genome_ver = get_arg("gv", "NA", \%args);
my $fixed_positive_color = get_arg("c", "", \%args);
my $fixed_positive_color_by_list_str = get_arg("cs", "", \%args);
my $fixed_negative_color = get_arg("mc", "", \%args);
my $zero_color_intensity = get_arg("zeroc", "", \%args);
my $min_color_intensity = get_arg("minc", "-10", \%args);
my $max_color_intensity = get_arg("maxc", "", \%args);
my $display_mode = get_arg("m", "By Types", \%args);
my $location_display_mode = get_arg("l", "Color gradient", \%args);
my $location_height = get_arg("lh", 10, \%args);
my $use_fixed_type_ordering = get_arg("fixed_order", 0, \%args);
my $display_track_information = get_arg("dt", 0, \%args);
my $add_feature_name_to_id = get_arg("af", 0, \%args);
my $add_counter_to_id = get_arg("ai", "", \%args);
my $feature_vector = get_arg("v", 0, \%args);
my $chromosome_starts = get_arg("chr_starts", "", \%args);
my $chromosome_ends = get_arg("chr_ends", "", \%args);
my $footer_pixels   = get_arg("fp", 5, \%args);

if (length($chromosome_starts) > 0) { $chromosome_starts = "ChromosomesStarts=\"$chromosome_starts\""; }
if (length($chromosome_ends) > 0) { $chromosome_ends = "ChromosomesEnds=\"$chromosome_ends\""; }

my $type = $feature_vector == 1 ? "ChromosomeFeatureVectorTrack" : "ChromosomeFeatureTrack";

if (length($desc) > 0) 
{ 
   $desc = &remove_illegal_xml_chars($desc);
   $desc =~ s/\t/ /g;
}

print "<GeneXPressChromosomeTrack Type=\"$type\" Name=\"$name\" Organism=\"$org\" GenomeVersion=\"$genome_ver\" Description=\"$desc\" DisplayMode=\"$display_mode\" LocationDisplayMode=\"$location_display_mode\" LocationHeight=\"$location_height\" FooterPixels=\"$footer_pixels\" $chromosome_starts $chromosome_ends";

if ($use_fixed_type_ordering == 1)
{
    print "DisplayByTypesOrder=\"true\" ";
}

if ($display_track_information == 1)
{
    print "DisplayTrackInfo=\"true\" ";
}

my %fixed_positive_color_by_list;
if (length($fixed_positive_color_by_list_str) > 0)
{
    my @row = split(/\;/, $fixed_positive_color_by_list_str);
    for (my $i = 0; $i < @row; $i += 2)
    {
	$fixed_positive_color_by_list{$row[$i]} = $row[$i + 1];
    }
}

my @features;
my %feature_types2id;
my @feature_types;
my $feature_counter = 0;

my $curr_type;
while(<$file_ref>)
{
  chop;

  my @row = split(/\t/);
  $curr_type = &remove_illegal_xml_chars($row[4]);
  if (length($feature_types2id{$curr_type}) == 0)
  {
    $feature_types2id{$curr_type} = $feature_counter++;
    push(@feature_types, $curr_type);
  }

  push(@features, $_);
}

print "NumLocations=\"".($#features+1)."\" FeatureTypes=\"";
for (my $i = 0; $i < @feature_types; $i++)
{
  if ($i > 0) { print ";" } 

  print "$feature_types[$i]";
}
print "\" ";

for (my $i = 0; $i < @feature_types; $i++)
{
  print "ColorGradientMaxPositiveColor_$i=\"" . &get_positive_color($i, $feature_types[$i]) . "\" ";
  print "ColorGradientMinNegativeColor_$i=\"" . &get_negative_color($i) . "\" ";

  if (length(get_arg("zeroc$i", "", \%args)) > 0)
  {
      print "ColorGradientZeroValue_$i=\"" . get_arg("zeroc$i", "", \%args) . "\" ";
  }
  elsif (length($zero_color_intensity) > 0)
  {
      print "ColorGradientZeroValue_$i=\"$zero_color_intensity\" ";
  }

  if (length(get_arg("minc$i", "", \%args)) > 0)
  {
      print "ColorGradientMinNegativeValue_$i=\"" . get_arg("minc$i", "", \%args) . "\" ";
  }
  elsif (length($min_color_intensity) > 0)
  {
      print "ColorGradientMinNegativeValue_$i=\"$min_color_intensity\" ";
  }

  if (length(get_arg("maxc$i", "", \%args)) > 0)
  {
      print "ColorGradientMaxPositiveValue_$i=\"" . get_arg("maxc$i", "", \%args) . "\" ";
  }
  elsif (length($max_color_intensity) > 0)
  {
      print "ColorGradientMaxPositiveValue_$i=\"$max_color_intensity\" ";
  }
}

print ">\n";

my $counter = 1;
for (my $i = 0; $i < @features; $i++)
{
    my @row = split(/\t/, $features[$i], 6);
    $curr_type = &remove_illegal_xml_chars($row[4]);
    if ($add_feature_name_to_id == 1) { $row[1] .= " ($feature_types2id{$curr_type})"; }
    if ($add_counter_to_id == 1) { $row[1] .= " (ID$counter)"; $counter++;}

    print "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$feature_types2id{$curr_type}\t$row[5]\n";
}

print "</GeneXPressChromosomeTrack>\n";

sub get_positive_color ()
{
  my ($color_index, $feature_name) = @_;

  my $res = $fixed_positive_color;

  if (length($res) == 0)
  {
      $res = get_arg("c$color_index", "", \%args);

      if (length($res) == 0 and length($fixed_positive_color_by_list_str) > 0)
      {
	  foreach my $str (keys %fixed_positive_color_by_list)
	  {
	      if ($feature_name =~ /$str/)
	      {
		  $res = $fixed_positive_color_by_list{$str};
		  last;
	      }
	  }
      }

      if (length($res) == 0)
      {
	  if ($color_index % 18 == 0)  { $res = "255,0,0,1"; }
	  if ($color_index % 18 == 1)  { $res = "0,0,255,1"; }
	  if ($color_index % 18 == 2)  { $res = "255,128,128,1"; }
	  if ($color_index % 18 == 3)  { $res = "128,0,64,1"; }
	  if ($color_index % 18 == 4)  { $res = "128,255,128,1"; }
	  if ($color_index % 18 == 5)  { $res = "50,1,131,1"; }
	  if ($color_index % 18 == 6)  { $res = "0,0,0,1"; }
	  if ($color_index % 18 == 7)  { $res = "0,128,255,1"; }
	  if ($color_index % 18 == 8)  { $res = "255,128,192,1"; }
	  if ($color_index % 18 == 9)  { $res = "128,128,192,1"; }
	  if ($color_index % 18 == 10) { $res = "255,128,64,1"; }
	  if ($color_index % 18 == 11) { $res = "255,255,128,1"; }
	  if ($color_index % 18 == 12) { $res = "128,0,255,1"; }
	  if ($color_index % 18 == 13) { $res = "128,128,0,1"; }
	  if ($color_index % 18 == 14) { $res = "0,89,45,1"; }
	  if ($color_index % 18 == 15) { $res = "0,255,128,1"; }
	  if ($color_index % 18 == 16) { $res = "0,255,0,1"; }
	  if ($color_index % 18 == 17) { $res = "128,255,255,1"; }
      }
  }

  return $res;
}

sub get_negative_color ()
{
    my ($color_index) = @_;
  
    my $res = $fixed_negative_color;

    if (length($res) == 0)
    {
	$res = get_arg("mc$color_index", "", \%args);

	if (length($res) == 0)
	{
	    $res = "0,255,0,1";
	}
    }

    return $res;
}

sub remove_illegal_xml_chars
{
  my $str = $_[0];
  $str =~ s/\&/&amp;/g;
  $str =~ s/\"/&quot;/g;
  $str =~ s/\'/&apos;/g;
  $str =~ s/\</&lt;/g;
  $str =~ s/\>/&gt;/g;
  $str =~ s/;/,/g;
  return $str;
}

__DATA__

/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/tab2feature_gxt.pl <file> 

    Creates a feature gxt file from a tab file

    -n <name>:          Name of the chromosome track (default: Features )
    -o <str>:           Organism
    -gv <str>:          Genome Version (default: NA)
    -d <desc>:          Track description (Tabs will be converted to spaces, default: empty)
 
    -c <str>:           Fix the max color for all features to be str (e.g., '255,0,0,1')
    -c<num>:            Fix the max color for feature <num> to be str (e.g., '255,0,0,1')
    -cs <list>:         Fix the max color for features that contain words from the list to colors
                        specified in the list. The format of list is: <str1;color1;...>
                        for specifying that features whose name contain str1 are set to color1
    -mc <str>:          Fix the min color for all features to be str (e.g., '255,0,0,1')
    -mc<num> <str>:     Fix the min color for feature <num> to be str (e.g., '255,0,0,1')

    -zeroc <num>:       Fix the zero color intensity for all features to be num
    -minc <num>:        Fix the min color intensity for all features to be num
    -maxc <num>:        Fix the max color intensity for all features to be num
    -minc<num1> <num2>: Fix the min color intensity for feature <num1> to be <num2>
    -maxc<num1> <num2>: Fix the max color intensity for feature <num1> to be <num2>
    -zeroc<num1> <num2>:Fix the zero color intensity for feature <num1> to be <num2>

    -m <str>:           Display mode (Full/Packed/Dense/By Types) (default: By Types)
    -l <str>:           Location Display mode (Color gradient/Filled box/Unfilled box/Filled oval/Unfilled oval/Directed filled box/Directed unfilled box) (default: Color gradient)
    -lh <num>:          Location height (default: 10)

    -fixed_order:       Use a fixed ordering on the feature types

    -dt:                Display track information
    -fp <num>:          Footer pixels (spacing below the track).

    -af:                Add the id of the feature to the name/id of the feature
    -ai:                Add a counter id to each instance

    -v:                 Feature vector

    -chr_starts <str>   Chromosome start locations (optional format: <chr_name>;<start>...)
    -chr_ends <str>     Chromosome end locations (optional format: <chr_name>;<end>...)

