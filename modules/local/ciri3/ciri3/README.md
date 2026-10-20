# CIRI3 notes

CIRI3 v3.0.1 is supplied by the Bioconda `ciri3` package and invoked through its
`ciri3` wrapper. The module sets `JAVA_TOOL_OPTIONS` because this wrapper does not
honour `JAVA_OPTS`.

CIRI3 does not expose a version option, so the version topic records the pinned
Bioconda package version, `3.0.1`. The package is licensed GPL-2.0-only.

The upstream v3.0.1 README documents detection output
as 13 tab-separated columns. Columns 2-4 are chromosome, 1-based inclusive
start, and end; column 5 is BSJ read count; column 11 is strand; and column 12
contains supporting BSJ read IDs separated by commas. `-UC` takes a tab-separated
list of chromosome, 1-based start, and end (not a BED file).

In `-UC` mode, the documented matrix sidecars are
`${output}.BSJ_Matrix` and `${output}.FSJ_Matrix` (for this module,
`${prefix}.ciri3.txt.BSJ_Matrix` and `${prefix}.ciri3.txt.FSJ_Matrix`). The upstream
README identifies each circRNA by `chr:start|end`; it does not document matrix
header/column semantics in the `-UC` section, so consumers should verify those
before assuming a sample-matrix layout.
