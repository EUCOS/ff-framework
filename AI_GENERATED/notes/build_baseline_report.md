# Build Baseline Report

This report records the first local build and smoke-test baseline for `EUCOS/ff-framework`.

All work described here preserves the project rule:

- original human-written framework code remains untouched;
- no files under `src/`, `include/`, `cmd/`, `data/`, `Makefile`, `lib/`, `bin/`, or `obj/` are intentionally modified;
- AI-generated notes and future helper material live under `AI_GENERATED/`.

## Environment

- Platform: Ubuntu 24.04.1 LTS under WSL2.
- Build tools installed through Ubuntu packages.
- `build-essential` is installed.
- `git` is installed.
- Working directory: `~/EUCOS/ff-framework`.

## Library Build

Command:

`make all`

Result: succeeded.

Confirmed output:

`lib/libff.a`

The library exists and is approximately 32 MB.

Observed baseline warnings:

- repeated overloaded-virtual warnings involving `State` and `MatterState`;
- one deprecated `bind2nd` warning in `src/Model.cpp`.

These are baseline warnings from the current source tree. They should not be treated as AI-induced regressions.

## D4_Gauge_Searcher Command Build

Command:

`make bin/D4_Gauge_Searcher`

Result: failed.

Observed failure:

`cmd/D4_Gauge_Searcher/main.cpp:7:10: fatal error: FF_Systematic_Gauge_Builder.hh: No such file or directory`

Interpretation:

`cmd/D4_Gauge_Searcher/main.cpp` appears to use an old header name:

`#include "FF_Systematic_Gauge_Builder.hh"`

The modern include tree appears to use names such as:

`#include "SystematicGaugeBuilder.h"`

This is an old command-level compatibility issue, not a failed library build.

Important: do not edit the original `cmd/D4_Gauge_Searcher/main.cpp` unless the research group explicitly decides to manually modernize it later.

## Framework_Tester Command Build

Command:

`make bin/Framework_Tester`

Result: succeeded.

Confirmed output:

`bin/Framework_Tester`

The executable exists and is approximately 4.4 MB.

## Framework_Tester Smoke Run

Command:

`./bin/Framework_Tester`

Result: ran to completion.

Observed output includes:

- `Begin program ... forgive me.`
- `Linearly independent.`
- `k_ij good.`
- `A 3 1 | A 3 1 | A 3 1 | D 4 1 | E 8 1 |`
- `Total matter representations: 34`
- `U(1)'s: 1`
- `ST SUSYs: 1`
- `You screwed up the matter rep equivalence classes.`

Interpretation:

`Framework_Tester` reaches its final matter-representation equivalence check.

The final message is hard-coded in `cmd/Framework_Tester/main.cpp` when `Equal_Matter_Rep_Classes(...)` returns false. Therefore this is a baseline smoke-test result, not an AI-induced regression.

This does not block gauge-scan optimization work, because the immediate project goal is gauge-group generation and scan performance, not matter-representation equivalence repair.

## Current Build Baseline Summary

Confirmed:

- library build works under WSL;
- `lib/libff.a` is produced;
- a modern command, `Framework_Tester`, builds;
- `Framework_Tester` runs to completion;
- the old `D4_Gauge_Searcher` command does not currently build because of an outdated `FF_*` include.

## Next Recommended Step

Create an additive AI-generated modern D4 gauge-scan wrapper or benchmark harness under:

`AI_GENERATED/benchmarks/`

The wrapper should:

- use modern headers;
- link against `lib/libff.a`;
- avoid modifying original files;
- initially target a tiny D=4 layer-1 order-2 gauge scan or even a dry-run/argument-parsing benchmark;
- compare output against baseline counts and gauge-group summaries once a baseline is established.

Original human-written code should remain untouched.
