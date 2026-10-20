# psirc runtime

These modules use Bioconda psirc 1.0.0 build 2, including its native `kallisto_b`
fork. The executable path is passed explicitly to `psirc`; reconstruction also
receives the complete BSJ output directory, including its internal mapping state.
Transcriptome preparation supplies the exon numbering required by the packaged
helper, including for GTFs without the optional `exon_number` attribute. Exons
must have a `transcript_id`; `gene_name` is optional.

The index command exits 255 even after successful indexing. `PSIRC_INDEX` accepts
that status only with the `Finished generating` marker and both nonempty indexes.

psirc has no version option, so the version topic records the pinned package
version. The kallisto fork's version is queried separately. Its non-commercial
license restriction is recorded in the module metadata and pipeline usage docs.

The synthetic reconstruction test creates paired-end reads with shared read IDs
and Illumina-style space-separated mate annotations, as required by psirc's read
extraction. It generates real BSJ state before testing full-length reconstruction.
