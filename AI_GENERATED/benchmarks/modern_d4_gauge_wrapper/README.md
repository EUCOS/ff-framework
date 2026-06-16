# Modern D4 Gauge Wrapper

This directory contains an AI-generated wrapper for the D=4 layer-1 gauge scan.

Purpose:

- reproduce the intent of the old `cmd/D4_Gauge_Searcher` command using modern framework header names;
- keep the original human-written command file untouched;
- provide a small benchmark/smoke-test entry point for future profiling work.

The original `cmd/D4_Gauge_Searcher/main.cpp` currently includes old `FF_*` headers and does not build against the modern include tree. This wrapper uses modern headers such as:

- `ModelBuilder.h`
- `SystematicGaugeBuilder.h`
- `GSOCoefficientMatrix.h`

## Build

From the repository root, after running `make all`:

```bash
g++ -Wall -pedantic -g -pg \
  -Iinclude \
  AI_GENERATED/benchmarks/modern_d4_gauge_wrapper/main.cpp \
  lib/libff.a \
  -o AI_GENERATED/benchmarks/modern_d4_gauge_wrapper/D4_Gauge_Searcher_Modern
./AI_GENERATED/benchmarks/modern_d4_gauge_wrapper/D4_Gauge_Searcher_Modern 2

## Run small test

```bash
./AI_GENERATED/benchmarks/modern_d4_gauge_wrapper/D4_Gauge_Searcher_Modern 2
```text id="c1dd98"
