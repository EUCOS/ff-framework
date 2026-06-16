# Benchmark Harnesses

Future benchmark scripts, drivers, or harnesses can live here when the research group wants additive profiling or comparison tooling.

Benchmark material in this directory should not modify original source files. It should link against or execute the existing framework as an external target, capture baseline output, and compare original output to AI-generated prototype output.

Recommended benchmark behavior:

- keep inputs small and deterministic at first
- record exact build flags and commands
- compare final counts and gauge-group summaries
- preserve baseline output separately from prototype output
- treat elapsed time as a performance metric, not a correctness metric

Any prototype optimization should remain separate from the human-written framework code until reviewed and manually integrated by the research group.

