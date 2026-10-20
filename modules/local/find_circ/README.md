# find_circ runtime

Bioconda find_circ 1.2 still requires Python 2.7. The latest compatible packaged
Bowtie2 is 2.5.4; newer Bowtie2 builds require a Zlib ABI incompatible with the
legacy pysam stack. Keep conda-forge before bioconda in channel priority.

The Python 2.7 dependency prevents native ARM64 builds. These legacy modules
remain AMD64-only; the pipeline's independently packaged Bowtie2 modules are not
subject to this version hold.
