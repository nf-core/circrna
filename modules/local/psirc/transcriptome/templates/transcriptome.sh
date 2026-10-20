#!/usr/bin/env bash
set -euo pipefail

# The packaged helper requires exon_number, an optional GTF attribute. It uses
# these values as unique keys, then sorts exons by genomic position itself.
perl - '${gtf}' > psirc.numbered.gtf <<'PERL'
use strict;
use warnings;
my %exons;
while (<>) {
    if (/^#/ || /^\\s*\$/) {
        print;
        next;
    }
    chomp;
    my @fields = split /\\t/, \$_, -1;
    if (@fields == 9 && \$fields[2] eq 'exon') {
        my (\$transcript) = \$fields[8] =~ /(?:^|;\\s*)transcript_id\\s+"([^"]+)"/;
        die "Exon is missing transcript_id\\n" unless defined \$transcript;
        my \$number = ++\$exons{\$transcript};
        \$fields[8] =~ s/(?:^|;\\s*)exon_number\\s+[^;]+;?/;/g;
        \$fields[8] =~ s/\\s*\$//;
        \$fields[8] .= ';' unless \$fields[8] =~ /;\$/;
        \$fields[8] .= qq( exon_number "\$number";);
    }
    print join("\\t", @fields), "\\n";
}
PERL

create_custom_transcriptome_fa.pl '${fasta}' psirc.numbered.gtf '${prefix}.${suffix}'
