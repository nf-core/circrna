# CIRI3 bundled JAR provenance

`bin/CIRI3_Java_1.8.0.jar` is the CIRI3 v3.0.1 Java 8 build from the upstream
release source, pinned to commit
[`73108c4478a7e691cb93360921b65203e0701516`](https://github.com/gyjames/CIRI3/tree/73108c4478a7e691cb93360921b65203e0701516).

- Upstream project/release: <https://github.com/gyjames/CIRI3/releases/tag/v3.0.1>
- Source path: `CIRI3_Java_1.8.0.jar`
- SHA-256: `fd2075523efde219fffc3bf48c0e12bd3b3fab5093a4c81a41e7273772b1718d`
- Upstream license: GPL-2.0; the accompanying upstream `LICENSE.txt` is retained
  alongside this notice.

The JAR is run on Java 17. Its upstream v3.0.1 README documents detection output
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
