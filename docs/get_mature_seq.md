# get_mature_seq.sh

`get_mature_seq.sh` operates by attempting to overlap circRNAs with the remaining biotypes in the filtered GTF file. If there are no overlapping features, the circRNA is considered intronic (ciRNA). If the circRNA perfectly overlaps exon boundaries, it is considered a circRNA. In the situation where a circRNA overlaps features but does not perfectly overlap exon boundaries, 2 scenarios are tested:

1. Within 200nt of an exon boundary: attempt to fit as circRNA with the largest (highest exon count) transcript.
2. Falls outside of 200nt of an exon boundary: treated as an EI-circRNA.

The script returns circRNAs in BED12 format, for sequence extraction using bedtools.
