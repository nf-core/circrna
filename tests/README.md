# Test configuration

Module and subworkflow tests must not depend on `conf/modules.config`. They must
remain runnable when that file is empty, even though Nextflow normally loads it
through the pipeline's `nextflow.config`.

Follow the nf-core pipeline template and scdownstream selector conventions in
`conf/modules.config`: use bare process names or aliases for general settings,
and `.*:SUBWORKFLOW:PROCESS` where the calling context matters. Do not add the
pipeline's root workflow name to every selector.

- Declare required tool arguments and output prefixes in the component's
  `tests/nextflow.config`, and load it explicitly with `config "./nextflow.config"`.
- Configure setup processes there too. A subworkflow test may explicitly include
  its dependencies' test configs instead of duplicating their settings.
- A test-specific config that replaces the suite config must include the common
  component test config when it needs those options.
- Assert emitted outputs where possible. If a test inspects an intermediate that
  is not emitted, configure its publishing explicitly in the test config.

`tests/nextflow.config` supplies shared test-data locations and test resources.
Pipeline integration tests continue to use `conf/modules.config` and may override
its options through their own test configs.
