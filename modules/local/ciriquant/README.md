# CIRIquant runtime

All three CIRIquant modules use the custom image built from this directory's
`Dockerfile`. It installs the [`nictru/CIRIquant` fork at
`e4916ca7b3370cef54d76ca162be64792d8c1b16`](https://github.com/nictru/CIRIquant/tree/e4916ca7b3370cef54d76ca162be64792d8c1b16),
plus version-pinned primary runtime dependencies.

This is an intentional exception to the Bioconda `ciriquant=1.1.3` package:
that package does not support this pipeline's single-end quantification path
and has runtime failures. Consequently, these modules do not declare Conda
environments or Conda locks. Architecture support must be verified by building
and running the image; the presence of architecture-independent Python sources
alone does not establish ARM64 support.

The AMD64 image is published as `docker.io/nicotru/ciriquant:45acbdd4f4c8b517`
and pinned by digest in the modules. It was smoke-tested with the container
entrypoint bypassed, as required by Nextflow. ARM64 is not currently published:
its StringTie 3.0.3 build requires HTSlib <1.23, conflicting with SAMtools 1.24's
HTSlib >=1.24 requirement. The [native ARM64 build log](https://wave.seqera.io/view/builds/bd-4d0d7774c9420691_1)
records this solver failure; AMD64 resolves the same pinned dependencies.

The equivalent frozen AMD64 SIF is published at
`oras://docker.io/nicotru/ciriquant:35228e054dcad0b0`.

The Dockerfile adapts the fork's two removed `decideTestsDGE` calls to edgeR's
current `decideTests` generic. The positive DE regression exercises this path
alongside the existing missing-output failure guard.
