use strict;

sub fileFriendlyName
{
   my ($name) = @_;

   $name =~ s/\//_/g;
   $name =~ s/\\/_/g;
   $name =~ s/"/_/g;
   $name =~ s/'/_/g;
   $name =~ s/\(/_/g;
   $name =~ s/\)/_/g;
   $name =~ s/\[/_/g;
   $name =~ s/\]/_/g;
   $name =~ s/\{/_/g;
   $name =~ s/\}/_/g;
   $name =~ s/ +/_/g;

   $name =~ s/_[_]+/_/g;

   return $name;
}

sub readFile
{
   my ($file) = @_;
   my @lines;
   while(my $line = <$file>)
   {
      chomp($line);
      push(@lines, $line);
   }
   return \@lines;
}

sub readFileName
{
   my ($file_name) = @_;
   my $file  = &openFile($file_name, "r");
   my $lines = &readFile($file);
   close($file);
   return $lines;
}

##------------------------------------------------------------------------
## \*FILE openFile ($string file_name=undef, $string purpose="r")
##------------------------------------------------------------------------
sub openFile
{
  my $file_name = shift;
  my $purpose   = shift;
  $file_name    = defined($file_name) ? $file_name : '-';
  $purpose      = defined($purpose) ? $purpose : 'r';
  my $file      = undef;
  if($file_name eq '-')
  {
    # Read-only files
    if($purpose =~ /^ *r/i)
    {
      $file = \*STDIN;
    }
    # Overwrite the file
    elsif($purpose =~ /^ *w/i or $purpose =~ /^ *a/i)
    {
      $file = \*STDOUT;
    }
  }
  else
  {
    # Read-only files
    if($purpose =~ /^ *r/i or not(defined($purpose)))
    {
      if($file_name =~ /\.gz$/ or $file_name =~ /\.Z$/)
      {
        $file_name = "zcat $file_name |";
      }
      else
      {
        $file_name = "<$file_name";
      }
    }
    # Overwrite
    elsif($purpose =~ /^ *w/i)
    {
      $file_name = ">$file_name";
    }
    # Append
    elsif($purpose =~ /^ *a/i)
    {
      $file_name = ">>$file_name";
    }

    if(not(open($file, $file_name)))
    {
       $file = undef;
    }
  }
  return $file;
}

sub printFilep
{
   my ($file_in, $file_out) = @_;
   $file_in  = defined($file_in)  ? $file_in  : \*STDIN;
   $file_out = defined($file_out) ? $file_out : \*STDOUT;

   while(my $line = <$file_in>)
   {
      print $file_out $line;
   }
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub GetNumColumns
{
   my ($file_name) = @_;

   my $rows = `head -1 $file_name | /nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/transpose.pl | wc -l`;

   $rows =~ /^ *([^ ]+)/;

   return $1;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub GetNumRows
{
   my ($file_name) = @_;

   my $rows = `wc -l $file_name`;

   $rows =~ /^ *([^ ]+)/;

   return $1;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub GetAllDirs
{
   return &getAllDirs(@_);
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getAllDirs
{
  my @dirs = @_;
  my @subdirs;
  foreach my $dir (@dirs)
  {
     my $dirp;
     if(opendir($dirp, $dir))
     {
        my @files = readdir($dirp);
        foreach my $file (@files)
        {
          my $path = $dir . '/' . $file;
          if ((-d $path) and (not($file =~ /^\.\.$/)) and (not($file =~ /^\.$/)))
          {
             push(@subdirs, $path);
          }
        }
        closedir($dirp);
     }
  }
  return @subdirs;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getAllFiles
{
  my @dirs = @_;
  my @paths;
  foreach my $dir (@dirs)
  {
    my $dirp;
    if(opendir($dirp,$dir))
    {
      my @files = readdir($dirp);
      foreach my $file (@files)
      {
        push(@paths, $dir . '/' . $file);
      }
      closedir($dirp);
    }
  }
  return @paths;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getAllFilesRecursively # ($follow_links,@dirs)
{
  my $follow_links = shift;
  my @dirs = @_;
  my @files = ();
  my @allEntries;
  my @entries;
  my $entry;
  my $subEntry;

  foreach $entry (@dirs)
  {
    if(-d $entry and ($follow_links or not(-l $entry)))
    {
      if(opendir(DIR,$entry))
      {
        @allEntries = readdir(DIR);
        @entries = ();
        while(@allEntries)
        {
          $subEntry = shift @allEntries;
          if($subEntry ne '..' and $subEntry ne '.')
          {
            if($entry =~ /\/$/)
            {
              push(@entries, $entry . $subEntry);
            }
            else
            {
              push(@entries, $entry . '/' . $subEntry);
            }
          }
        }
        closedir(DIR);
        push(@files, &getAllFilesRecursively($follow_links,@entries));
      }
      else
      {
        print STDERR "Could not open directory -->$entry<-- skipping.\n";
      }
    }
    elsif(-f $entry)
    {
      push(@files,$entry);
    }
  }
  return @files;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getAllCodeRecursively # ($follow_links,@dirs)
{
  my $follow_links = shift;
  my @dirs = @_;
  my @files = ();
  my @allEntries;
  my @entries;
  my $entry;
  my $subEntry;

  foreach $entry (@dirs)
  {
    if(-d $entry and ($follow_links or not(-l $entry)))
    {
      if(opendir(DIR,$entry))
      {
        @allEntries = readdir(DIR);
        # print STDERR "[[$allEntries[0] $allEntries[1]]]\n";
        @entries = ();
        while(@allEntries)
        {
          $subEntry = shift @allEntries;
          if($subEntry ne '..' and $subEntry ne '.')
          {
            if($entry =~ /\/$/)
            {
              push(@entries, $entry . $subEntry);
            }
            else
            {
              push(@entries, $entry . '/' . $subEntry);
            }
          }
        }
        closedir(DIR);
        push(@files, &getAllCodeRecursively($follow_links,@entries));
      }
      else
      {
        print STDERR "Could not open directory -->$entry<-- skipping.\n";
      }
    }
    elsif(-f $entry and isCodeFile($entry))
    {
      push(@files,$entry);
    }
  }
  return @files;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getPathPrefix
{
  my $path = shift @_;

  $path =~ s/[^\/]*[\/]*$//;
  return $path;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getPathSuffix
{
  my $path = shift @_;

  while($path =~ /\/$/)
  {
    chop($path);
  }
  if($path =~ /([^\/]+)$/)
  {
    $path = $1;
  }
  return $path;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub remPathExt
{
  my $path = shift @_;

  $path =~ s/\.[^\.]*$//;

  return $path;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getPathExt
{
  my $path = shift @_;

  my $ext = '';
  if($path =~ /(\.[^\.]*)$/)
  {
    $ext = $1;
  }

  return $ext;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub expandPath
{
  my $file = shift @_;
  my $home = "$ENV{HOME}";

  # print STDERR "'$file' --> ";
  $file =~ s/^ +//;
  $file =~ s/ +$//;
  $file =~ s/~/$home/ge;
  $file =~ s/\$HOME/$home/ge;
  $file =~ s/\$\(HOME\)/$home/ge;
  $file =~ s/\$\{HOME\}/$home/ge;
  # print STDERR "'$file'\n";

  return $file;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub resolvePath # ($file,$pwd)
{
  my $file = &expandPath(shift);
  my $pwd  = shift;
  $pwd =~ s/[\/]+ *$//;

  $file = &isRelative($file) ? ($pwd . '/' . $file) : $file;

  return $file;
}

##------------------------------------------------------------------------
## Returns the depth of the file from the root directory.
##------------------------------------------------------------------------
sub getPathDepth
{
  my $file = &expandPath(shift);
  # Remove any leading ./ that indicate "current" directory.
  while($file =~ s/^\.\///) {}

  my $depth=0;
  for(; $file =~ /\S/; $file = &getPathPrefix($file))
    { $depth++; }

  return $depth;
}


##------------------------------------------------------------------------
## Returns 1 if the file is an absolute path, 0 otherwise.
##------------------------------------------------------------------------
sub isAbsolute # ($file)
{
  my $file = &expandPath(shift);

  # If the first symbol is /, it's absolute.
  return ($file =~ /^ *\//);
}

##------------------------------------------------------------------------
## Returns 1 if the file is a relative path, 0 otherwise.
##------------------------------------------------------------------------
sub isRelative # ($file)
{
  my $file = shift;
  return not(&isAbsolute($file));
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub dos2unix
{
  my $file = shift @_;

  $file =~ s/\\/\//g;
  $file =~ s/[Cc]:/\/cygdrive\/c/;

  return $file;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub isCodeFile
{
  my $file = shift;

  if(not(-f $file))
  {
    return 0;
  }

  if($file =~ /\.c$/i)            { return 1; }
  if($file =~ /\.cc$/i)           { return 1; }
  if($file =~ /\.cpp$/i)          { return 1; }
  if($file =~ /\.h$/i)            { return 1; }
  if($file =~ /\.hh$/i)           { return 1; }
  if($file =~ /\.htm$/i)          { return 1; }
  if($file =~ /\.html$/i)         { return 1; }
  if($file =~ /\.hpp$/i)          { return 1; }
  if($file =~ /\.java$/i)         { return 1; }
  if($file =~ /\.mak$/i)          { return 1; }
  if($file =~ /\.map$/i)          { return 1; }
  if($file =~ /\.pl$/i)           { return 1; }
  if($file =~ /\.pm$/i)           { return 1; }
  if($file =~ /\.py$/i)           { return 1; }
  if($file =~ /\.m$/i)            { return 1; }
  if($file =~ /\.s$/i)            { return 1; }
  if($file =~ /\.y$/i)            { return 1; }
  if($file =~ /Makefile/)         { return 1; }
  if($file =~ /Makefile\.common/) { return 1; }

  return 0;
}


##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub isObjectFile
{
  my $file = shift;

  if(not(-f $file))
  {
    return 0;
  }

  if($file =~ /\.o$/i) { return 1; }
  if($file =~ /\.so$/i) { return 1; }

  return 0;
}


##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub protectFromShell # ($file)
{
  my $file = shift;
  $file =~ s/([ ~"')(&])/\\$1/g;
  return $file;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub deprotectFromShell # ($file)
{
  my $file = shift;
  $file =~ s/\\([^\\])/$1/g;
  return $file;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub getFileText # ($file_name)
{
  my $file_name = shift;

  if (defined($file_name))
  {
    my $file      = &openFile($file_name,'r');
    my @text      = <$file>;
    my $text      = join("",@text);
    return $text;
  }
  else
  {
    return "";
  }
}

##------------------------------------------------------------------------
## $int numRows ($string file_name, $int blanks=1)
##------------------------------------------------------------------------
sub numRows
{
   my ($file_name, $blanks) = @_;

   $blanks = defined($blanks) ? $blanks : 1;

   my $num_rows = 0;

   if(-f $file_name)
   {
      my $wc = $blanks ? (`wc $file_name`) : (`grep -v '^[ ][ ]*\$' $file_name | wc`);

      my @wc = split(/ +/, $wc);

      $num_rows = $wc[1];
   }
   return $num_rows;
}

sub pad
{
   my ($list, $filler) = @_;

   my @padded;

   for(my $i = 0; $i < scalar(@{$filler}); $i++)
   {
      push(@padded, defined($$list[$i]) ? $$list[$i] : $$filler[$i]);
   }

   return \@padded;
}

##------------------------------------------------------------------------
## $string fill ($int max_cols, \@list filling=[''], $int pad=0,
##               $string delim="\t", $string line=$_)
##------------------------------------------------------------------------
sub fill
{
   my ($max_cols, $filling, $pad, $delim, $line) = @_;
   $pad     = defined($pad)     ? $pad     :    0;
   $delim   = defined($delim)   ? $delim   : "\t";
   $line    = defined($line)    ? $line    :   $_;

   my $chomped = 0;
   if($line =~ /\n$/)
   {
      chomp($line);
      $chomped = 1;
   }

   my @tuple    = split($delim, $line);
   my $num_cols = scalar(@tuple);
   my $num_fill = defined($filling) ? scalar(@{$filling}) : 0;

   if(not($pad))
   {
      for(my $i = 0; $i < $num_cols; $i++)
      {
         if(length($tuple[$i]) == 0 and $num_fill > 0)
         {
            $tuple[$i] = $$filling[$i % $num_fill];
         }
      }
   }

   for(my $i = $num_cols; $i < $max_cols; $i++)
   {
      my $filler = $num_fill > 0 ? $$filling[$i % $num_fill] : '';
      push(@tuple, $filler);
   }

   my $filled = join($delim, @tuple);

   if($chomped)
   {
      $filled .= "\n";
   }

   return $filled;
}

##------------------------------------------------------------------------
## $int numTokens ($string delim="\t", $string line=$_)
##------------------------------------------------------------------------
sub numTokens
{
   my ($delim, $line) = @_;
   $delim = defined($delim) ? $delim : "\t";
   $line  = defined($line)  ? $line  : $_;
   return scalar(split($delim, $line));
}

##------------------------------------------------------------------------
## $int numCols ($string file='-', $string delim="\t", \$string eaten=undef)
##
## Takes a file name as an argument.  If the file is STDIN then it has to 
## eat a line to determine the columns.  This eaten line is saved in
## eaten if supplied.
##------------------------------------------------------------------------
sub numCols
{
   my ($file, $delim, $eaten) = @_;
   $file  = defined($file)  ? $file  : '-';
   $delim = defined($delim) ? $delim : "\t";

   my $filep;
   open($filep, $file) or die("Could not open file '$file'");
   my $num_cols = &numCols2($filep, $delim, $eaten);

   if($file ne '-')
   {
      close($filep);
      $$eaten = undef;
   }

   return $num_cols;
}

##-----------------------------------------------------------------------------
## $int numCols2 (\*FILE file=\*STDIN, $string delim="\t", \$string eaten=undef)
##
## Takes a file pointer as an argument.
##
## Eats a line of the file to figure out how many columns there are.
##-----------------------------------------------------------------------------
sub numCols2
{
   my ($file, $delim, $eaten) = @_;

   $file        = defined($file)  ? $file  : \*STDIN;
   $delim       = defined($delim) ? $delim : "\t";

   my $line     = <$file>;
   my $num_cols = &numTokens($delim, $line);

   if(defined($eaten))
   {
      $$eaten = $line;
   }

   return $num_cols;
}

##------------------------------------------------------------------------
## $int maxCols ($string file='-', $string delim="\t", \@list eaten=undef)
##
## Takes a file name as an argument.  If the file is STDIN then it has to 
## eat a line to determine the columns.  This eaten line is saved in
## eaten if supplied.
##------------------------------------------------------------------------
sub maxCols
{
   my ($delim, $file, $eaten) = @_;

   $file  = defined($file)  ? $file  : '-';

   $delim = defined($delim) ? $delim : "\t";

   my $filep;

   open($filep, $file) or die("Could not open file '$file'");

   my $max_cols = 0;

   while(my $line = <$filep>)
   {
      my $num_cols = &numTokens($delim, $line);

      if($num_cols > $max_cols)
      {
         $max_cols = $num_cols;
      }

      if($file eq '-' and defined($eaten))
      {
         push(@{$eaten}, $line);
      }
   }
   close($filep);

   return $max_cols;
}

##----------------------------------------------------------------------------
## \@list &readHeader(\*FILE file=\*STDIN, $int num_rows=1, $string delim="\t")
#
#  Read columns from the file.  Returns a list containing each column
#  where the elements in the column are seperated by delim.
##----------------------------------------------------------------------------
sub readHeader
{
   my ($file, $num_rows, $delim_in, $delim_out) = @_;
   $file      = defined($file)      ? $file      : \*STDIN;
   $num_rows  = defined($num_rows)  ? $num_rows  : 1;
   $delim_in  = defined($delim_in)  ? $delim_in  : "\t";
   $delim_out = defined($delim_out) ? $delim_out : $delim_in;

   my @column;

   for(my $i = 0; ($i < $num_rows) and not(eof($file)); $i++)
   {
      my $line  = <$file>;

      my @tuple = split($delim_in, $line);

      chomp($tuple[$#tuple]);

      for(my $j = 0; $j < scalar(@tuple); $j++)
      {
         $column[$j] .= ($i == 0) ? $tuple[$j] : ($delim_out . $tuple[$j]);
      }
   }

   return \@column;
}

sub duplicate
{
   my ($num_copies, @entries) = @_;

   $num_copies = defined($num_copies) ? $num_copies : 1;

   if(scalar(@entries) == 0)
   {
      push(@entries, '');
   }

   my @tuple;

   for(my $i = 0; $i < $num_copies; $i++)
   {
      push(@tuple, @entries);
   }

   return \@tuple;
}

##----------------------------------------------------------------------------
## \@list &getCols(\@list, $string ranges, \$int max_cols)
##----------------------------------------------------------------------------
sub getCols
{
   my ($list, $ranges, $max_cols) = @_;

   if(defined($list) and defined($ranges) and defined($max_cols))
   {
      my $num_cols  = scalar(@{$list});
      my $updated   = 0;
      if($num_cols > $$max_cols)
      {
         $$max_cols = $num_cols;
         if(defined($ranges))
         {
            my @cols = &parseRanges($ranges, $$max_cols);
            for(my $i  = 0; $i < scalar(@cols); $i++)
            {
               $cols[$i]--;
            }
            return \@cols;
         }
      }
   }
   return undef;
}


##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub fileExists
{
   my ($file) = @_;
   return (-f $file);
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub dirExists
{
   my ($dir) = @_;
   return (-d $dir);
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub linkExists
{
   my ($link) = @_;
   return (-l $link);
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub pathExists
{
   my ($path) = @_;
   return ((-f $path) or (-d $path) or (-l $path));
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub allFilesExist
{
   foreach my $file (@_)
   {
      if(not(&fileExists($file)))
      {
         return 0;
      }
   }
   return 1;
}

##------------------------------------------------------------------------
##
##------------------------------------------------------------------------
sub allPathsExist
{
   foreach my $path (@_)
   {
      if(not(&pathExists($path)))
      {
         return 0;
      }
   }
   return 1;
}

##------------------------------------------------------------------------
## \@list getHeader ($string file, $int num_lines=1, $string delim="\t"
##                   \$FILE* filep)
##
##    \$FILE* filep - If supplied gets a copy of the file pointer (and
##                    the file is *not* closed but kept open).
##------------------------------------------------------------------------
sub getHeader
{
   my ($file, $num_lines, $delim, $filep) = @_;
   $num_lines = defined($num_lines) ? $num_lines : 1;
   $delim     = defined($delim) ? $delim : "\t";

   my $fp;
   open($fp, $file) or die("Could not open file '$file'");

   my @header;
   for(my $i = 0; $i < $num_lines; $i++)
   {
      my $line = <$fp>;
      my @line = split($delim, $line);
      chomp($line[$#line]);

      for(my $j = 0; $j < scalar(@line); $j++)
      {
         $header[$j] .= length($header[$j]) > 0 ? " $line[$j]" : $line[$j];
      }
   }

   if(not(defined($filep)))
   {
      close($fp);
   }
   else
   {
      $$filep = $fp;
   }

   return \@header;
}

#------------------------------------------------------------------------
# @list parseRanges ($string ranges, $int last_col, $int offset=0)
#------------------------------------------------------------------------
sub parseRanges
{
   my ($ranges, $last_col, $offset) = @_;
   $offset = defined($offset) ? $offset : 0;

   my @segments                     = ();
   my @range                        = ();
   my @fields                       = ();
   my ($i,$beg,$end,$inc)           = (0,-1,-1,1);
   my $seg;
   @segments = split(",", $ranges);
   foreach $seg (@segments)
   {
      $beg = undef;
      $end = undef;
      $inc = undef;

      if($seg =~ /^([-]{0,1}\d+)[-:](\d+)[-:]([-\d]*)/)
      {
         $beg = $1;
         $inc = $2;
         $end = $3;
      }
      elsif($seg =~ /^([-]{0,1}\d+)[:-]([-\d]*)/)
      {
         $beg = $1;
         $end = $2;
         $inc = 1;
      }
      elsif($seg =~ /^([-\d]+)/)
      {
         $beg = $1;
         $end = $1;
         $inc = 1;
      }

      $beg = (defined($beg) and $beg =~ /\S/) ? $beg : 0;
      $beg = (defined($beg) and $beg < 0)     ? (defined($last_col) ? $last_col + $beg + 1 : undef) : $beg;
      $end = (defined($end) and $end =~ /\S/) ? $end : (defined($last_col) ? $last_col : undef);
      $end = (defined($end) and $end < 0)     ? (defined($last_col) ? $last_col + $end + 1 : undef) : $end;
      $inc = $inc =~ /\S/ ? $inc : 1;

      if(defined($beg) and defined($end) and defined($inc))
      {
         for($i = $beg; $i <= $end; $i += $inc)
         {
            push(@fields, $i);
         }
      }
   }

   if($offset != 0)
   {
      for(my $i = 0; $i < scalar(@fields); $i++)
      {
         $fields[$i] += $offset;
      }
   }

   return @fields;
}

##--------------------------------------------------------------------------------------
## $int evalRegex ($string value=$_, $string regex='/\S/',
##                    $string op='=~', $string hardop=undef, $int negate=0)
##
## Returns 1 if the regular expression evaluates to true.
##
## hardop - supplies a specific operation.  Valid are:
##             'gt'  for $value >  $regex
##             'gte' for $value >= $regex
##             'lt'  for $value <  $regex
##             'lte' for $value <= $regex
##--------------------------------------------------------------------------------------
sub evalRegex
{
   my ($value, $regex, $op, $hardop, $negate) = @_;
   $value  = defined($value)  ? $value  : $_;
   $regex  = defined($regex)  ? $regex  : '/\S/';
   $op     = defined($op)     ? $op     : '=~';
   $hardop = defined($hardop) ? $hardop : undef;
   $negate = defined($negate) ? $negate : 0;

   my $result = 0;

   if(not(defined($hardop)))
   {
      $result = $negate ? eval("not(\$value $op $regex);") :
                          eval("\$value $op $regex;");
      $result = (defined($result) and $result) ? 1 : 0;
   }
   elsif($hardop eq 'eq')
   {
      $result = $value eq $regex ? 1 : 0;
   }
   elsif($hardop eq 'ne')
   {
      $result = $value ne $regex ? 1 : 0;
   }
   elsif($hardop eq 'gt')
   {
      $result = $value gt $regex ? 1 : 0;
   }
   elsif($hardop eq 'lt')
   {
      $result = $value lt $regex ? 1 : 0;
   }
   elsif($hardop eq 'gte')
   {
      $result = (($value gt $regex) or ($value eq $regex)) ? 1 : 0;
   }
   elsif($hardop eq 'lte')
   {
      $result = (($value lt $regex) or ($value eq $regex)) ? 1 : 0;
   }
   elsif($hardop eq '=')
   {
      $result = $value == $regex ? 1 : 0;
   }
   elsif($hardop eq '==')
   {
      $result = $value == $regex ? 1 : 0;
   }
   elsif($hardop eq '!=')
   {
      $result = $value != $regex ? 1 : 0;
   }
   elsif($hardop eq '>')
   {
      $result = $value > $regex ? 1 : 0;
   }
   elsif($hardop eq '>=')
   {
      $result = $value >= $regex ? 1 : 0;
   }
   elsif($hardop eq '<')
   {
      $result = $value < $regex ? 1 : 0;
   }
   elsif($hardop eq '<=')
   {
      $result = $value <= $regex ? 1 : 0;
   }
   return $result;
}

##------------------------------------------------------------------------
## \%attrib parseArgs (\@list args, \@list flag_tuples, $int store_extra=0)
##
## \@args        - list of arguments (from @ARGV usually)
##
## \@flag_tuples - a list of triples where each triple gives:
##
##                        (FLAG,TYPE,DEFAULT,EATS)
##
##           for each option.  FLAG is the string (or regular expression)
##           that indicates use of the option (including the dash if 
##           dashes are being used!).  DEFAULT is the default value for
##           the option used when the option's flag is not found in the
##           argument list.  If EATS is defined then the next argument in
##           \@args is "eaten" to retrieve the value for the option; if
##           EATS is undef then it does *not* eat the next character
##           but fills in the string contained in EATS for the option
##           when its flag is encountered. TYPE can be any of:
##
##                'scalar' - the option is a scalar quantity
##                'list'   - the option is a list
##                'set'    - the option is a set
##                'map'    - the option is an associative array
##                'file'   - the option is the name of a file
##
##
##           The special FLAG '--file' reads a file from the arguments.
##
## $store_extra - If 1 stores unrecognized arguments in a list that can
##                be accessed with '--extra'.  Otherwise the function
##                dies if an unrecognized argument is encountered (default)
##
##------------------------------------------------------------------------
sub parseArgs
{
   my ($args, $flag_tuples, $store_extra) = @_;
   $store_extra = defined($store_extra) ? $store_extra : 0;

   my @args = defined($args)        ? @{$args}        : ();
   my @f    = defined($flag_tuples) ? @{$flag_tuples} : ();
   my %options;
   my %option_overwritten;
   my @extra;

   # Store the defaults in the options.
   for(my $i = 0; $i < scalar(@f); $i++)
   {
      my ($flag, $type, $default, $eats) = ($f[$i][0],$f[$i][1],$f[$i][2],$f[$i][3]);
      $options{$flag} = $default;
   }

   while(@args)
   {
      my $arg = shift @args;
      if($arg eq '--help')
      {
         @args    = ();
         $options{'--help'} = 1;
      }
      else
      {
         my $matched_arg = 0;
         for(my $i = 0; $i < scalar(@f); $i++)
         {
            my ($flag, $type, $default, $eats) = ($f[$i][0],$f[$i][1],$f[$i][2],$f[$i][3]);

            if(($flag =~ /^--file/ and ((-f $arg) or ($arg eq '-'))) or
               ($arg eq $flag))
            {
               $matched_arg = 1;
               my $value = $arg;
               if(defined($eats))
               {
                  $value = $eats;
               }
               elsif(not($flag =~ /^--file/))
               {
                  $value = shift @args;
               }

               if($type eq 'scalar')
               {
                  $options{$flag} = $value;
               }
               elsif($type eq 'list')
               {
                  if(not(defined($option_overwritten{$flag})))
                  {
                     my @empty_list;
                     $options{$flag} = \@empty_list;
                  }
                  my $list = $options{$flag};
                  push(@{$list}, $value);
               }
               elsif($type eq 'set')
               {
                  if(not(defined($option_overwritten{$flag})))
                  {
                     my %empty_set;
                     $options{$flag} = \%empty_set;
                  }
                  my $set = $options{$flag};
                  $$set{$value} = 1;
               }
               elsif($type eq 'map')
               {
                  if(not(defined($option_overwritten{$flag})))
                  {
                     my %empty_map;
                     $options{$flag} = \%empty_map;
                  }
                  my ($attrib,$val) = split('=',$value);
                  my $map = $options{$flag};
                  $$map{$attrib} = $val;
               }
               $option_overwritten{$flag} = 1;
            }
         }
         if(not($matched_arg))
         {
            if($store_extra)
            {
               push(@extra, $arg);
            }
            else
            {
               die("Bad argument '$arg' given");
            }
         }
      }
   }

   # If we picked up any extra arguments put them in the hash.
   $options{'--extra'} = \@extra;

   return \%options;
}

sub readKeyedTuples
{
   my ($file, $delim, $key_col, $headers) = @_;
   $file    = defined($file)    ? $file    : $_;
   $delim   = defined($delim)   ? $delim   : "\t";
   $key_col = defined($key_col) ? $key_col : 0;
   $headers = defined($headers) ? $headers : 0;

   my %set;

   my $filep   = &openFile($file);
   my $line_no = 0;
   while(<$filep>)
   {
      $line_no++;
      if($line_no > $headers)
      {
         my @x = split($delim);
         chomp($x[$#x]);
         my $key = splice(@x, $key_col, 1);
         $set{$key} = \@x;
      }
   }

   return \%set;
}

sub download
{
   my ($url, $dir, $pattern) = @_;
   $dir = defined($dir) ? $dir : '.';

   my $wget = "wget --retr-symlinks --passive-ftp --follow-ftp -Q 0 -N -nd -t 1 -P $dir";

   my $cmd  = defined($pattern) ? "$wget '$url/*' -A$pattern" :
                                  "$wget '$url'";

   `$cmd`;

   my $downloaded_file = undef;

   if($url =~ /\/([^\/]+) *$/)
   {
      $downloaded_file = $dir . '/' . $1;
   }

   return $downloaded_file;
}

sub downloadMultiple
{
   my ($url, $dir, $pattern) = @_;
   $dir = defined($dir) ? $dir : '.';

   my $wget = "wget --retr-symlinks --passive-ftp -Q 0 -N -nd -t 1 -P $dir";

   my $cmd  = "$wget $url";

   `$cmd`;

   my $downloaded_file = undef;

   if($url =~ /\/([^\/]+) *$/)
   {
      $downloaded_file = $dir . '/' . $1;
   }

   return $downloaded_file;
}

sub copyFile
{
   my ($source, $destination) = @_;

   my $s = &openFile($source, "r");

   my $d = &openFile($destination, "w");

   (defined($s) and defined($d)) or die("Could not copy '$source' to '$destination'");

   while(my $line = <$s>)
   {
      print $d $line;
   }

   close($s);

   close($d);
}

sub deleteFile
{
   my ($file) = @_;

   system("rm -f $file");
}

sub makeTmpFile
{
   my ($prefix, $file) = @_;

   $prefix = defined($prefix) ? $prefix : '';

   my $suffix = '.tmp';

   my $tmpname  = '';

   my $tmpfile  = undef;

   my $maxtries = 1000000;

   for(my $i = 0; not(defined($tmpfile)) and ($i < $maxtries); $i++)
   {
      my $time = time;

      my $rand = int(rand(1000000));

      $tmpname = $prefix . '_' . $time . '_' . $rand . $suffix;

      if(&pathExists($tmpname))
      {
         $tmpfile = undef;
      }

      else
      {
         if(defined($file))
         {
            &copyFile($file, $tmpname);
         }
         return $tmpname;
      }
   }

   die("Could not open a temporary file");

   return $tmpname;
}

sub findLine
{
   my ($filep, $patterns, $logic) = @_;

   $logic = defined($logic) ? $logic : 'and';

   my $result = undef;

   if(($logic =~ /and/i) or ($logic =~ /all/i))
   {
      $result = &findLineWhereAllMatch($filep, $patterns);
   }

   if(($logic =~ /or/i) or ($logic =~ /any/i))
   {
      $result = &findLineWhereAnyMatch($filep, $patterns);
   }

   return $result;
}

sub findLineWhereAllMatch
{
   my ($filep, $patterns) = @_;

   while(my $line = <$filep>)
   {
      my $match = 1;

      foreach my $pattern (@{$patterns})
      {
         if(not($line =~ /$pattern/))
         {
            $match = 0;
         }
      }

      if($match)
      {
         return $line;
      }
   }

   return undef;
}

sub findLineWhereAnyMatch
{
   my ($filep, $patterns) = @_;

   while(my $line = <$filep>)
   {
      my $match = 0;

      foreach my $pattern (@{$patterns})
      {
         if($line =~ /$pattern/)
         {
            $match = 1;
         }
      }

      if($match)
      {
         return $line;
      }
   }

   return undef;
}

sub printMatrix
{
   my ($matrix, $file, $delim_col, $delim_row) = @_;

   $file  = defined($file)  ? $file  : \*STDOUT;

   $delim_col = defined($delim_col) ? $delim_col : "\t";

   $delim_row = defined($delim_row) ? $delim_row : "\n";

   my $n = scalar(@{$matrix});

   for(my $i = 0; $i < $n; $i++)
   {
      print $file join($delim_col, @{$$matrix[$i]}), $delim_row;
   }
}

sub readIds
{
   my ($file, $column, $delim) = @_;

   $column = defined($column) ? $column : 0;

   $delim  = defined($delim)  ? $delim  : "\t";

   my $filep = &openFile($file);

   my $num_keys = 0;

   my %seen;

   my @ids;

   my @rows;

   my $row = 0;

   while(my $line = <$filep>)
   {
      my @tuple = split($delim, $line);

      if($column <= $#tuple)
      {
         my $key = $tuple[$column];

         if(not(exists($seen{$key})))
         {
            $ids[$num_keys]  = $key;

            $rows[$num_keys] = [$row];

            $seen{$key} = $num_keys;

            $num_keys++;
         }
         else
         {
            my $first_row = $seen{$key};

            push(@{$rows[$first_row]}, $row);
         }
      }

      $row++;
   }

   close($filep);

   return (\@ids, \@rows);
}

sub readDataMatrix
{
   my ($file, $key_col, $delim, $max_cols) = @_;

   $file    = defined($file) ? $file : $_;

   $key_col = defined($key_col) ? $key_col : 0;

   $delim   = defined($delim) ? $delim : "\t";

   my $filep = &openFile($file);

   my @data;

   if(defined($max_cols))
   {
      $$max_cols = 0;
   }

   while(my $line = <$filep>)
   {
      my @tuple = split($delim, $line);

      chomp($tuple[$#tuple]);

      my $key   = splice(@tuple, $key_col, 1);

      unshift(@tuple, $key);

      push(@data, \@tuple);

      if(defined($max_cols) and scalar(@tuple) > $$max_cols)
      {
         $$max_cols = scalar(@tuple);
      }
   }

   close($filep);

   return \@data;
}

sub getFilesFromList
{
   my ($args) = @_;

   my @files;

   foreach my $arg (@{$args})
   {
      if((-f $arg) or (-l $arg))
      {
         push(@files, $arg);
      }
   }

   return \@files;
}

sub getColumn
{
   my ($file, $col, $delim, $blank) = @_;

   $file  = defined($file) ? $file : $_;

   $col   = defined($col)  ? $col  : 0;

   $delim = defined($delim) ? $delim : "\t";

   my @list;

   my $filep = &openFile($file);

   while(my $line = <$filep>)
   {
      my @x = split($delim, $line);

      chomp($x[$#x]);

      my $val = $col < scalar(@x) ? $x[$col] :
                (defined($blank) ? $blank : undef);

      if(defined($val))
      {
         push(@list, $val);
      }
   }

   close($filep);

   return \@list;
}


##--------------------------------------------------------------------------------------
## ExtractFileLinesBetweenRegexps
##--------------------------------------------------------------------------------------
sub ExtractFileLinesBetweenRegexps
{
  my ($in_file_ref, $out_file_ref, $regexp_from, $regexp_to, $from_inclusive, $to_inclusive) = @_;

  if ( $regexp_from eq $regexp_to ) {
    die "(/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl) ExtractFileLinesBetweenRegexps() - 'from' and 'to' regexps are identical. Give me a break, use grep.\n";
  }

  my $reached_from = 0;
  if ( $regexp_from eq "" ) {
    $reached_from = 1;
  }

  while (<$in_file_ref>) {
    chomp;
    
    if ( $reached_from == 0 ) {
      if ( $_ =~ /$regexp_from/ ) {
	$reached_from = 1;
	next unless ( $from_inclusive );
      }
    }
    
    if ( $reached_from == 1 ) {
      if ( $regexp_to ne "" and $_ =~ /$regexp_to/ ) {
	print $out_file_ref "$_\n" if ( $to_inclusive );
	last; # break from while loop, and end.
      }
      else {
	print $out_file_ref "$_\n";
      }
    }
  }
}


##--------------------------------------------------------------------------------------
## ExtractAllFileLinesNotBetweenRegexps
##--------------------------------------------------------------------------------------
sub ExtractAllFileLinesNotBetweenRegexps
{
  my ($in_file_ref, $out_file_ref, $regexp_from, $regexp_to, $from_inclusive, $to_inclusive) = @_;

  if ( $regexp_from eq $regexp_to ) {
    die "(/nfs/data3/CIRCEST/pipeline_sponging/modules/local/pita/templates/lib/libfile.pl) ExtractAllFileLinesNotBetweenRegexps() - 'from' and 'to' regexps are identical. Give me a break, use grep.\n";
  }

  my $print_flag = 1;
  if ( $regexp_from eq "" ) {
    $print_flag = 0;
  }

  while (<$in_file_ref>) {
    chomp;

    if ( $print_flag == 1 ) {
      if ( $_ =~ /$regexp_from/ ) {
	$print_flag = 0;
	next unless ( $from_inclusive );
      }
      print $out_file_ref "$_\n";
    }

    if ( $print_flag == 0 ) {
      last if ( $regexp_to eq "" ); # break from while loop, and end.

      if ( $_ =~ /$regexp_to/ ) {
	$print_flag = 1;
	print $out_file_ref "$_\n" if ( $to_inclusive );
      }
    }
  }
}


##--------------------------------------------------------------------------------------
## ReplaceStringInFile
##--------------------------------------------------------------------------------------
sub ReplaceStringInFile
{
  my ($in_file_ref, $out_file_ref, $to_replace, $replace_with) = @_;

  while (<$in_file_ref>) {
    chomp;
    $_ =~ s/$to_replace/$replace_with/g;
    print $out_file_ref "$_\n";
  }
}

1
