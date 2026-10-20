# Nanopore runtime

The current AMD64 environment uses SAMtools 1.24. A native ARM64 build cannot
resolve this stack: circtools pulls StringTie, whose ARM64 builds require HTSlib
<1.23, while SAMtools 1.24 requires HTSlib >=1.24.

The [ARM64 build report](https://wave.seqera.io/view/builds/bd-40f295e234e21811_1)
records the dependency conflict. This module is therefore AMD64-only; the other
circtools modules have separately resolved environments.
