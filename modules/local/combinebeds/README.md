# Plotting compatibility

The `filter` and `reads` modules retain pandas 2.3.3 alongside UpSetPlot 0.9.0.
UpSetPlot's default-style handling uses chained in-place assignments that fail
under pandas 3's Copy-on-Write semantics, producing invalid `NaN` plot colors.
Both modules keep the current Polars and Matplotlib releases; other environments
can use pandas 3 independently.

Overlap plots and their embedded MultiQC image sections remain required outputs
when there is more than one tool or sample to compare.
