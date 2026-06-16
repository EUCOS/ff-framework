# Code Map

This note summarizes the current repository structure and the apparent model-building pipeline, based on read-only inspection.

## Repository Structure

- `Makefile`: builds the static library and command binaries.
- `src/`: framework implementation files.
- `include/`: framework headers.
- `cmd/`: command-line programs, each with a `main.cpp`.
- `AI_GENERATED/`: AI-generated notes and future additive helper material.

Some command files use older `FF_*` include/type names, while the current `include/` tree uses names such as `ModelBuilder.h`, `SystematicBuilder.h`, and `SystematicGaugeBuilder.h`. That compatibility question should be checked before treating command build failures as framework regressions.

## Command Entry Points

The `cmd/` tree contains searchers, model openers/builders, statistics generators, and comparison tools. For systematic gauge scans, the apparent command-level entry point is `cmd/D4_Gauge_Searcher/main.cpp`, which sets `Large_ST_Dimensions = 4`, loads the S vector and default `k_ij`, constructs a gauge-specific systematic builder, and calls `Perform_Search()`.

`cmd/Framework_Tester/main.cpp` uses the current header names and is a useful smoke-test candidate.

## SystematicGaugeBuilder

`SystematicGaugeBuilder` is the gauge-specialized scan path. Its `Perform_Search()` method writes front matter, builds `k_ij` extensions if needed, creates the extension basis-vector slots, and starts recursive layer construction.

Its main scan work is split across:

- `Build_Extensions()`: recursive layer driver.
- `Build_Basis_Vectors()`: builds gauge chunks and checks chunk combinations for modular invariance.
- `Build_Models()`: loops over `k_ij` extensions, builds gauge-group-only models, and inserts `LEEFT` summaries into the uniqueness set.

## SystematicBuilder

`SystematicBuilder` is the generic systematic scan base. It stores the initial model, extension orders, common basis alphas, counters, output streams, `k_ij` extensions, and `Unique_LEEFTs_`.

It also provides:

- `Build_k_ij_Extensions()`
- `Form_k_ij_Layer_Extensions()`
- generic basis-vector extension building
- generic full model construction

## Basis-Vector Generation

Basis vectors are represented by `BasisVector`. Systematic generation flows through `MIBVGenerator`, which builds left-moving chunks, compact right-moving chunks, observable chunks, and hidden chunks.

Important helper generators include:

- `LMGenerator`
- `RealChunkGenerator`
- `ComplexChunkGenerator`

Gauge scans use `MIBVGenerator::Build_Gauge_Chunks()`, which fixes the left-moving gauge chunk and generates the right-moving chunks needed for gauge-model enumeration.

## Modular-Invariance Checks

Full model modular invariance is checked by `ModularInvarianceChecker`.

During systematic scans, chunk-level pruning uses `ChunkConsistencyChecker`, especially:

- `Check_Modular_Invariance()`
- `Check_D10_Modular_Invariance()`
- `Check_Simultaneous_Periodic_Modes()`

These checks sit inside high-volume nested loops and are likely performance-critical.

## k_ij, GSO, and GGSO Matrix Construction

The half-matrix data is stored in `GSOCoefficientMatrix`.

Related builders are:

- `GSOCoefficientMatrixGenerator`: enumerates possible extension rows.
- `GSOCoefficientMatrixBuilder`: converts and completes a half GSO matrix, computes diagonal/off-diagonal elements, and checks consistency.

The code commonly refers to `k_ij` and GSO coefficient matrix data. Any GGSO terminology should be mapped carefully to the existing `GSOCoefficientMatrix` and `GSOProjector` implementation before changing behavior.

## ModelBuilder

`ModelBuilder` is the main construction class. It loads basis vectors and `k_ij` rows, checks consistency, builds basis alphas, builds alpha sectors, builds the fermion-mode map, completes the GSO matrix, and constructs gauge groups and matter content.

Gauge-only scans call `Build_Gauge_Group_Model()`, which avoids the full matter-state path and builds:

- basis alphas
- alpha sectors
- fermion-mode map
- complete `k_ij`
- gauge groups
- U(1) factors

Full model scans call `Build_Model()`, which additionally builds SUSY states and matter states.

## StateBuilder

`StateBuilder` builds physical states for a sector. It delegates massless state generation to left- and right-moving builders, combines candidates, and applies GSO projection.

Important paths:

- `Build_Boson_States()`
- `Build_Fermion_States()`
- `Build_SUSY_States()`

## GSOProjector

`GSOProjector` applies GSO projections to candidate states through:

- `GSOP_Boson()`
- `GSOP_Fermion()`

These functions loop over common basis alphas, state numerators, coefficients, and fermion-mode map entries. They are likely relevant in full model scans and gauge-group construction.

## GaugeGroupIdentifier

`GaugeGroupIdentifier` takes positive roots and identifies gauge group class, rank, Kac-Moody level, and simple roots.

Important work includes:

- class identification across ADE and BCGF families
- simple-root search
- connectivity checks
- simple-root ordering for representation conventions

This can be expensive when many positive roots are generated.

## LEEFT Uniqueness Filtering

Systematic scans store unique low-energy effective field theory summaries in `std::set<LEEFT>`.

Uniqueness is controlled by `LEEFT::operator<`, `LEEFT::operator==`, and matter-representation equivalence checks through `MatterRepEquivalenceTester`. This can become a cost center because every insertion may trigger ordered-set comparisons and equivalence testing.

