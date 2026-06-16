# Correctness Baselines

Run small, deterministic baselines before any optimization or prototype comparison. The purpose is to catch behavior changes early while scans are still small enough to inspect.

## Build-Only Baseline

First confirm the library build independently:

```sh
make clean
make all
```

Then check command builds separately. Some command files appear to use older `FF_*` headers and type names, so command build failures should be triaged as compatibility issues unless the specific command is already known to build in the current environment.

## Framework_Tester Smoke Test

`cmd/Framework_Tester/main.cpp` uses current header names and exercises several important paths:

- NAHE loading
- basis-vector loading
- linear-independence checks
- `k_ij` consistency checks
- modular-invariance checks
- full model construction
- particle-content display
- `LEEFT` construction and comparison
- matter-representation equivalence checks

Use it as a smoke test before comparing scan outputs.

## D4_Gauge_Searcher Order 2

The intended small gauge-scan baseline is:

```sh
D4_Gauge_Searcher 2
```

For this baseline, compare:

- final total unique models
- final total consistent models
- final total consistent basis vectors
- final total basis vectors tested
- gauge-group summaries in model output
- generated model-extension basis vectors and `k_ij` rows where deterministic

Do not compare elapsed time as a correctness condition.

## D4_Gauge_Searcher Order 3

If applicable and tractable in the current environment, also run:

```sh
D4_Gauge_Searcher 3
```

Use the same comparison criteria as the order-2 gauge scan. This gives coverage for odd-order behavior without relying on the no-supersymmetry path.

## No Odd-Order N=0 Behavior

Avoid using odd-order no-supersymmetry behavior as an early correctness baseline. It is too easy to mix intended physics constraints, scan assumptions, and implementation details. Establish lower-risk order-2 and applicable order-3 gauge baselines first.

## What To Compare

Prefer stable semantic outputs:

- final count lines
- gauge-group names and U(1) counts
- `k_ij`/GSO rows written for unique models
- basis-vector extension data
- `LEEFT` uniqueness counts

Avoid treating these as correctness signals:

- elapsed time
- order of debug-only messages unless intentionally deterministic
- profiler output
- file modification timestamps

## Baseline Storage

If saved baselines are needed later, store them only in an approved additive location such as `AI_GENERATED/benchmarks/` or another explicitly approved generated-output directory. Do not modify original source or data files to store baseline results.

