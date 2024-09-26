#!/usr/bin/perl
use File::Basename;
use strict;

require "${moduleDir}/templates/lib/load_args.pl";
require "${moduleDir}/templates/lib/format_number.pl";

my %args = load_args(\\@ARGV);

my \$prefix = 'circRNA';
my \$utr_fn = '${fasta}';
my \$upstream_fn = get_arg("upstream", "", \\%args);
my \$mir_fn = '${mature}';
my \$flank_up = get_arg("flank_up", "", \\%args);
my \$flank_down = get_arg("flank_down", "", \\%args);
my \$lengths = '${params.pita_lengths}';
my \$gu = '${params.pita_wobbles}';
my \$mismatches = '${params.pita_mismatches}';
my \$loop = get_arg("loop", "", \\%args);
my \$gxp = get_arg("gxp", "", \\%args);

if (length(\$utr_fn) == 0)
{
   die ("Error: must specifiy a UTR file (-utr)");
}
if (length(\$mir_fn) == 0)
{
   die ("Error: must specifiy a miRs file (-mir)");
}
if (length(\$prefix) == 0)
{
   die ("Error: must specifiy a prefix for the output file(s) (-prefix)");
}

my \$pid = \$\$;
my \$gxp_header = "<?xml version=\\"1.0\\" encoding=\\"iso-8859-1\\"?>\\n<GeneXPress>\\n\\n<TSCRawData>\\nUID\\tNAME\\tGWEIGHT\\tE1\\nG2\\tG2\\t1\\t1\\n</TSCRawData>\\n\\n<GeneXPressAttributes>\\n<Attributes Id=\\"0\\">\\n<Attribute Name=\\"g_module\\" Id=\\"0\\" Value=\\"1 2\\">\\n</Attribute>\\n</Attributes>\\n</GeneXPressAttributes>\\n\\n<GeneXPressObjects>\\n<Objects Type=\\"Genes\\" URLPrefix=\\"genome-www4.stanford.edu/cgi-bin/SGD/locus.pl?locus=\\">\\n<Gene Id=\\"0\\" ORF=\\"G2\\" Desc=\\"G2\\">\\n<Attributes AttributesGroupId=\\"0\\" Type=\\"Full\\" Value=\\"1\\">\\n</Attributes>\\n</Gene>\\n</Objects>\\n<Objects Type=\\"Experiments\\">\\n<Experiment Id=\\"0\\" name=\\"E1\\">\\n</Experiment>\\n</Objects>\\n</GeneXPressObjects>\\n<TSCHierarchyClusterData NumClusters=\\"1\\">\\n<Root ClusterNum=\\"0\\" NumChildren=\\"0\\">\\n</Root>\\n</TSCHierarchyClusterData>\\n";

print STDERR "Calculating predictions...\\n";

my \$runcmd = "${moduleDir}/templates/lib/pita_run.pl -utr \$utr_fn -mir \$mir_fn -prefix \${prefix}_ " . 
	(length(\$upstream_fn) > 0 ? " -upstream \$upstream_fn" : "") .
	(length(\$flank_up) > 0 ? " -flank_up \\"\$flank_up\\"" : "") .
	(length(\$flank_down) > 0 ? " -flank_down \\"\$flank_down\\"" : "") .
        (length(\$lengths) > 0 ? " -l \\"\$lengths\\"" : "") .
	(length(\$gu) > 0 ?  " -gu \\"\$gu\\"" : "") .
	(length(\$loop) > 0 ? " -loop \\"\$loop\\"" : "") .
	(length(\$mismatches) > 0 ? " -m \\"\$mismatches\\"" : "");

#print STDERR "Running: \$runcmd\\n";

system (\$runcmd);


system "rm \${prefix}_ext_utr.stab; mv \${prefix}_pita_results.tab tmp_\${pid}; head -n 1 tmp_\${pid} | cut -f 1-5,8,10- > \${prefix}_pita_results.tab; cat tmp_\${pid} | ${moduleDir}/templates/lib/body.pl 2 -1 | tr -d '\\r' | cut -f 1-8,10- | ${moduleDir}/templates/lib/merge_columns.pl -1 4 -2 5 -d ':' | ${moduleDir}/templates/lib/merge_columns.pl -1 4 -2 5 -d ':' | sed 's/Seed:Mismatchs:G:U/Seed/g' | sort -k 13n >> \${prefix}_pita_results.tab; rm tmp_\${pid};";

system "cat \${prefix}_pita_results.tab | ${moduleDir}/templates/lib/body.pl 2 -1 | cut -f 1,2,13 | ${moduleDir}/templates/lib/modify_column.pl -c 2 -m '\\"-1\\"' | ${moduleDir}/templates/lib/average_rows.pl -k 0,1 -losoe -n | ${moduleDir}/templates/lib/cut.pl -f 2,3,1,4 | ${moduleDir}/templates/lib/modify_column.pl -c 3 -m '\\"-1\\"' | ${moduleDir}/templates/lib/modify_column.pl -c 3 -p 2 | sort -k 4n | ${moduleDir}/templates/lib/cap.pl \\"RefSeq,microRNA,Sites,Score\\" > \${prefix}_pita_results_targets.tab";

print STDERR "Done.\\n";

if (\$gxp == 1)
{
   print STDERR "Creating gxp file...\\n";
   
   open (GXP_FILE, ">\${prefix}_pita_results.gxp");
   print GXP_FILE \$gxp_header;
   close GXP_FILE;

   system "${moduleDir}/templates/lib/fasta2stab.pl \$utr_fn | ${moduleDir}/templates/lib/stab2length.pl | ${moduleDir}/templates/lib/lin.pl | ${moduleDir}/templates/lib/add_column.pl -s 0 | ${moduleDir}/templates/lib/cut.pl -f 2,1,4,3 | ${moduleDir}/templates/lib/add_column.pl -s 0 | ${moduleDir}/templates/lib/add_column.pl -s 1 | ${moduleDir}/templates/lib/tab2feature_gxt.pl -n 'UTR Sequences' >> \${prefix}_pita_results.gxp";
   system "${moduleDir}/templates/lib/body.pl 2 -1 \${prefix}_pita_results.tab | tr -d '\\r' | ${moduleDir}/templates/lib/cut.pl -f 1-4,2,14 | ${moduleDir}/templates/lib/uniquify.pl -c 1 | ${moduleDir}/templates/lib/tab2feature_gxt.pl -n 'PITA Predictions' -c '0,0,255,1'  -lh 50 -l 'Filled box' -minc '0' -maxc '25' >> \${prefix}_pita_results.gxp";


   open (GXP_FILE, ">>\${prefix}_pita_results.gxp");
   print GXP_FILE "<GeneXPressTable Type=\\"ChromosomeTrack\\" Name=\\"UTR Sequences Track\\" TrackNames=\\"UTR Sequences\\">\\n</GeneXPressTable>\\n<GeneXPressTable Type=\\"ChromosomeTrack\\" Name=\\"PITA Predictions Track\\" TrackNames=\\"PITA Predictions\\">\\n</GeneXPressTable>\\n<TableDisplay TableDataModel=\\"PITA Predictions Track\\">\\n</TableDisplay>\\n<TableDisplay TableDataModel=\\"UTR Sequences Track\\">\\n</TableDisplay>\\n<GeneXPressClusterLists>\\n</GeneXPressClusterLists>\\n</GeneXPress>\\n";
   close GXP_FILE;

   print STDERR "Done.\\n";
}

sub removeIllegalXMLChars
{
   my \$str = \$_[0];
   
   my \$res_str = "";
   for (my \$i = 0; \$i < length(\$str); \$i++) {
      my \$char = substr(\$str, \$i, 1);
      if ((ord(\$char) >= 32 and ord(\$char) <= 126) or ord(\$char) == 10 or ord(\$char) == 9) {
	 \$res_str .= \$char;
	 }
   }
   
   \$res_str =~ s/\\&/&amp;/g;
   \$res_str =~ s/\\"/&quot;/g;
   \$res_str =~ s/\\'/&apos;/g;
   \$res_str =~ s/\\</&lt;/g;
   \$res_str =~ s/\\>/&gt;/g;
   
   return \$res_str;
}

__DATA__
syntax: pita_prediction.pl [OPTIONS]

Execute the PITA algorithm for identifying and scoring microRNA target sites.

options:
    -utr <filename>:      fasta file containing the UTRs to be scanned
    -mir <filename>:      fasta file containing the microRNA sequences
    -upstream <filename>: fasta file containing the upstream sequence for each UTR. The IDs
                          in should match the IDs found int the UTR file. If less 200 bp are
                          given (or if no file is given), it is padded with Poly-A.

    -flank_up <bp>
    -flank_down <bp>:     Flank requirement in basepairs (default: zero for both)
    
    -ddG_context <bp>:    Number of bases upstream and downstream for target site that are
                          taken into account when folding the UTR (default: 70)


    -prefix <string>:     Add the string as a prefix to the output files (pita_results.tab and ext_utr.stab)
    -gxp:        Produce a gxp (Genomica project file) output file.

    Seed matching parameters:

    -l <num1-num2>:       Search for seed lengths of num1,...,num2 to the MicroRNA (default: 6-8)

    -gu <nums>:           Lengths for which G:U wobbles are allowed and number of allowed wobbles.
                          Format of nums: <length;num G:U>,<length;num G:U>,... (default: 6;0,7;1,8;1)

    -m <nums>:            Lengths for which mismatches are allowed and number of allowed mismatches
                          Format of nums: <length;num mismatches>,<length;num mismatches>,...
                          (default: 6;0,7;0,8;1)

    -loop <nums>:         Lengths for which a single loop in either the target or the microrna is allowed
                          Format of nums: <length>,<length>,... (default: none)


