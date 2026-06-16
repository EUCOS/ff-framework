# Profiling Plan

This plan avoids modifying original source files. Profiling should use compiler flags, external tools, and additive harnesses only.

## Goals

- Identify true runtime hotspots before optimization.
- Separate build compatibility issues from performance findings.
- Keep all profiling support outside the human-written framework code.
- Preserve baseline output so optimized or prototype paths can be compared safely.

## Build Variants Without Source Edits

Use command-line Makefile variable overrides rather than changing `Makefile`.

Examples:

```sh
make clean
make exec FLAGS="-Wall -pedantic -g -O2 -pg"
make release FLAGS="-Wall -pedantic -O3 -g -fno-omit-frame-pointer"
```

The exact command should be selected after confirming which command programs currently compile with the available headers.

## External Profiling Options

On Linux or WSL:

- `gprof` with `-pg` builds.
- `perf record` and `perf report` with frame pointers.
- `valgrind --tool=callgrind` for call counts and instruction-level attribution.
- `callgrind_annotate` or KCachegrind/QCachegrind for report inspection.

On Windows:

- Visual Studio Performance Profiler.
- Windows Performance Recorder and Windows Performance Analyzer.
- Sampling profilers that can consume symbols from a debug build.

The preferred first pass is sampling, because systematic scans may spend time across many nested calls and containers. Instrumentation can follow only after sampling identifies narrow targets.

## Candidate Hotspots To Validate

Expected hotspots from read-only inspection:

- chunk generation in `MIBVGenerator`, `LMGenerator`, `RealChunkGenerator`, and `ComplexChunkGenerator`
- chunk modular-invariance checks in `ChunkConsistencyChecker`
- recursive extension construction in `SystematicBuilder` and `SystematicGaugeBuilder`
- per-candidate `ModelBuilder` reconstruction
- alpha-sector enumeration in `AlphaBuilder`
- state construction and projection in `StateBuilder` and `GSOProjector`
- gauge-root grouping and identification in `GaugeGroupIdentifier`
- `LEEFT` uniqueness insertion and matter-representation equivalence checks

These are hypotheses until measured.

## Additive Benchmark Harness

If profiling needs a controlled driver, add a new benchmark harness under `AI_GENERATED/benchmarks/` or another explicitly approved additive location. Do not edit existing `src/`, `include/`, `cmd/`, `Makefile`, `data/`, `lib/`, `bin/`, or `obj/` files.

A future harness can:

- include current public framework headers
- link against `lib/libff.a`
- construct small `ModelBuilder` initial conditions
- run one-layer gauge scans for small orders
- write results under an AI-generated output directory
- compare final counts and gauge-group summaries against saved baseline output

The harness should be disposable and clearly marked as AI-generated research support.

## Reporting

Profiling reports should include:

- exact branch and commit
- compiler and flags
- command or harness invocation
- input orders and scan type
- wall-clock time as context only
- final model counts and summaries for correctness
- profiler top functions and call stacks

Elapsed time should never be used as the only correctness signal.

