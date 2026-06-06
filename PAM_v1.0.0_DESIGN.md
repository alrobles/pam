# pam v1.0.0 — Architecture & Implementation Design

## Vision

Transform the simple R-only `pam` package into a high-performance,
CRAN-ready tool for constructing presence-absence matrices in **both
geographic and environmental space**, with a C++17 engine via Rcpp.

The problem is deceptively simple — convert a stack of species range
polygons into a binary sites × species matrix — but at scale (5000+
species, millions of cells, dual-space projection) it demands compiled
code, parallelism, sparse storage, and careful algorithmic choices.

------------------------------------------------------------------------

## 1. Current State (v0.2.0)

    get_geo_pam(shp, res)
      └─ terra::rast() → terra::rasterize() → terra::as.data.frame() → tidyr::spread()

**Bottlenecks:** -
[`terra::rasterize()`](https://rspatial.github.io/terra/reference/rasterize.html)
is single-threaded C++ (terra’s own, not GDAL) - Rasterizes ALL species
in one call via a categorical field — fast for \<100 species but chokes
on 5000+ with overlapping ranges - Returns dense tibble (not sparse) —
5000 species × 65000 cells = 325M entries - No environmental space
support - No coverage fraction or area-weighting option - No parallel
processing

------------------------------------------------------------------------

## 2. Proposed Architecture

### 2.1 Dual-Space PAM Framework

                        ┌──────────────────────────┐
                        │   Input: SpatVector       │
                        │   (species range polygons │
                        │    + "sciname" column)     │
                        └─────────┬────────────────┘
                                  │
                  ┌───────────────┼───────────────┐
                  ▼                               ▼
        ┌──────────────────┐            ┌──────────────────┐
        │  Geographic PAM  │            │ Environmental PAM │
        │  (G-PAM)         │            │ (E-PAM)           │
        │                  │            │                   │
        │  Grid: lat/lon   │            │  Grid: PCA axes   │
        │  res in degrees  │            │  of bioclim vars  │
        │  or km           │            │  (via xclim)      │
        └──────────────────┘            └──────────────────┘
                  │                               │
                  ▼                               ▼
        ┌──────────────────┐            ┌──────────────────┐
        │  C++ Engine:     │            │  C++ Engine:      │
        │  pip_rasterize() │            │  env_classify()   │
        │  (point-in-poly  │            │  (bin occurrence  │
        │   + R-tree index)│            │   points into     │
        │                  │            │   E-space grid)   │
        └────────┬─────────┘            └────────┬─────────┘
                 │                                │
                 └────────────┬───────────────────┘
                              ▼
                    ┌──────────────────┐
                    │  Sparse PAM      │
                    │  (dgCMatrix or   │
                    │   dense tibble)  │
                    └──────────────────┘

### 2.2 Module Breakdown

| Module | Language | Purpose |
|----|----|----|
| `src/pam_rasterize.cpp` | C++17 | Point-in-polygon rasterization with OpenMP |
| `src/pam_rtree.h` | C++17 | Bounding-box spatial index (header-only) |
| `src/pam_env.cpp` | C++17 | Environmental space binning + PCA projection |
| `src/pam_sparse.cpp` | C++17 | Sparse triplet → dgCMatrix conversion |
| `src/pam_omp.h` | C++17 | OpenMP thread safety (from xclim pattern) |
| `R/geo_pam.R` | R | `geo_pam()` — geographic PAM constructor |
| `R/env_pam.R` | R | `env_pam()` — environmental PAM constructor |
| `R/pam_class.R` | R | S3 class `pam` with print/summary/plot methods |
| `R/optimal_resolution.R` | R | Resolution optimizer for E-space |

------------------------------------------------------------------------

## 3. C++ Engine Design

### 3.1 Geographic PAM: Point-in-Polygon Approach

**Why not use terra::rasterize?** - terra rasterizes all species via a
single categorical field — breaks with overlapping ranges -
Single-threaded - Returns full raster, not sparse cell IDs

**Our approach: cell-centroid PIP (point-in-polygon)**

For macroecological PAMs at ≥0.5° resolution, a centroid-based PIP test
is standard practice (Arita et al. 2008, 2012; letsR package). The
algorithm:

    1. Generate grid of cell centroids for the target extent + resolution
    2. For each species polygon:
       a. Compute bounding box → filter candidate centroids (spatial index)
       b. For each candidate centroid: ray-casting PIP test
       c. Record (cell_id, species_id) pairs as sparse triplets
    3. Assemble sparse matrix from triplets

**C++ pseudocode:**

``` cpp
// Core: parallel PIP rasterization
// Input:  cell_x[], cell_y[]  — centroid coordinates (n_cells)
//         poly_coords[]       — packed polygon vertices per species
//         species_id[]        — species index per polygon
// Output: sparse triplets (cell_id, species_id)

struct Triplet { int cell; int species; };

std::vector<Triplet> pip_rasterize(
    const double* cx, const double* cy, int n_cells,
    const PolygonPack& polys, int n_species, int ncores)
{
    std::vector<std::vector<Triplet>> thread_local_results(ncores);

    #pragma omp parallel for schedule(dynamic) num_threads(ncores)
    for (int s = 0; s < n_species; ++s) {
        int tid = omp_get_thread_num();
        auto& bbox = polys.bbox(s);
        for (int c = 0; c < n_cells; ++c) {
            if (!bbox.contains(cx[c], cy[c])) continue;
            if (point_in_polygon(cx[c], cy[c], polys.ring(s))) {
                thread_local_results[tid].push_back({c, s});
            }
        }
    }
    // merge thread-local vectors
    return merge(thread_local_results);
}
```

**PIP algorithm: winding number** - More robust than ray-casting for
edge cases (vertices on boundary) - O(V) per point where V = polygon
vertices - Standard implementation from Computational Geometry in C
(O’Rourke)

**Spatial index strategy:** - Simple approach: per-species bounding-box
filter on cell centroids - For very large grids: lightweight R-tree over
cell centroids using a header-only implementation (no Boost dependency
for CRAN)

### 3.2 Environmental PAM

**Algorithm:**

    1. Extract bioclimatic values at cell centroids (or occurrence points)
       using terra::extract() or xclim rasters
    2. Run PCA on the environmental matrix (base R prcomp)
    3. Define E-space grid:
       - Compute range of PC1, PC2 (optionally PC3)
       - Divide into n_bins per axis → n_bins² cells in E-space
    4. For each species: classify occurrence points into E-space bins
    5. Record (e_cell_id, species_id) presence pairs
    6. Assemble E-PAM as sparse matrix

**Optimal resolution in E-space:**

Following Broennimann et al. (2012) and ecospat conventions: - Default:
kernel-density smoothing to avoid empty-cell artifacts - For discrete
PAM: use Sturges’ rule or Scott’s rule on PC scores
`n_bins = ceil(1 + log2(n_points))` per axis - Provide
`auto_resolution()` function that optimizes bin count to maximize
information content while avoiding sparsity:
`optimal_bins = argmin_k { KL_div(PAM_k, PAM_smooth) + λ·sparsity(k) }`

### 3.3 Sparse Storage

The PAM is inherently sparse (most species absent from most cells). For
5000 species × 65000 cells at 1° resolution, occupancy ≈ 5-15%.

**Strategy:** - Internal: sparse triplet vector `(i, j, x)` in C++ -
Export to R as `Matrix::dgCMatrix` (compressed sparse column) -
Convenience method `as_tibble()` for backward compatibility with
v0.2.0 - Parquet export via `arrow::write_parquet()` for HPC pipelines

### 3.4 Zero-Copy Bridge (from xclim pattern)

Following the `bioclim_xt()` pattern in xclim: - Input: `REAL()` pointer
aliasing from R `NumericMatrix` — no copy - Output: pre-allocate R
matrix, write directly into `REAL()` pointer - OpenMP: respect
`OMP_THREAD_LIMIT` via `pam_safe_threads()` helper (copy of
`xclim_omp.h` pattern)

------------------------------------------------------------------------

## 4. API Design

### 4.1 Core Functions

``` r

# Geographic PAM (replaces get_geo_pam)
geo_pam(
  x,                    # SpatVector with 'sciname' column
  res       = 1.0,      # grid resolution in CRS units
  extent    = NULL,      # optional: terra::ext() object
  crs       = NULL,      # optional: target CRS (default: from x)
  method    = "pip",     # "pip" (C++) or "terra" (legacy)
  sparse    = TRUE,      # return dgCMatrix (TRUE) or tibble (FALSE)
  ncores    = 1L,        # OpenMP threads
  ...
)

# Environmental PAM (NEW)
env_pam(
  x,                    # SpatVector with 'sciname' column
  env,                  # SpatRaster of environmental layers (e.g. bioclim)
  n_axes    = 2L,       # number of PCA axes to use
  n_bins    = "auto",   # bins per axis ("auto", or integer)
  method    = "pca",    # "pca" or "raw" (no PCA, use layers directly)
  sparse    = TRUE,
  ncores    = 1L,
  ...
)

# Unified constructor
pam(
  x,
  space     = "geo",    # "geo" or "env"
  ...                   # passed to geo_pam() or env_pam()
)
```

### 4.2 S3 Class `pam`

``` r
# Structure
pam_obj <- list(
  matrix   = <dgCMatrix or tibble>,   # the PAM itself
  space    = "geo" | "env",           # which space
  coords   = <tibble: x, y | PC1, PC2>,  # site coordinates
  species  = <character vector>,       # species names (column order)
  res      = <numeric>,               # resolution used
  crs      = <character>,             # CRS (geo only)
  pca      = <prcomp object>,         # PCA rotation (env only)
  meta     = list(                    # metadata
    n_sites   = ...,
    n_species = ...,
    fill      = ...,  # fraction of 1s
    created   = Sys.time()
  )
)
class(pam_obj) <- "pam"

# Methods
print.pam(x)          # compact summary
summary.pam(x)        # richness stats, range-size distribution
plot.pam(x)           # richness map (geo) or E-space occupancy (env)
as_tibble.pam(x)      # convert to dense tibble (backward compat)
richness(x)           # per-site species richness (row sums)
range_size(x)         # per-species range size (column sums)
```

### 4.3 Backward Compatibility

``` r

# Deprecated but still works:
get_geo_pam(shp, res = 0.5, ...)
# → internally calls geo_pam(shp, res, sparse = FALSE) + tibble conversion
# → emits .Deprecated() warning pointing to geo_pam()
```

------------------------------------------------------------------------

## 5. Dependencies

### DESCRIPTION (v1.0.0)

``` yaml
Package: pam
Title: Presence-Absence Matrix Construction in Geographic and Environmental Space
Version: 1.0.0
Depends: R (>= 4.1.0)
Imports:
    Rcpp (>= 1.0.0),
    checkmate,
    methods,
    terra,
    Matrix,
    tibble
LinkingTo:
    Rcpp
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown,
    arrow,
    ggplot2,
    sf
SystemRequirements: C++17
```

**Rationale:** - `Rcpp` — C++ bridge (same as maxentcpp, nicher, sobol,
xclim) - `Matrix` — sparse dgCMatrix for efficient PAM storage - `terra`
— polygon input handling + environmental raster extraction - `checkmate`
— input validation (lab standard) - Drop: `dplyr`, `tidyr`, `readr` (no
longer needed with C++ engine) - `arrow` in Suggests only (parquet
export for HPC users)

------------------------------------------------------------------------

## 6. Implementation Roadmap

### Phase 1: C++ Core — Geographic PAM (v0.3.0)

**SSDLC Stage 1: R reference implementation**

    geo_pam_r(x, res) → pure base R, no Rcpp
      - Generate grid centroids with seq()
      - terra::relate(vect, vect, "contains") for PIP
      - Build dense matrix with which()

**SSDLC Stage 2: C++ implementation**

    geo_pam_cpp(cx, cy, poly_x, poly_y, poly_idx, n_species, ncores)
      - Winding-number PIP in C++17
      - OpenMP parallelization over species
      - Returns integer triplets (i, j, 1)

**SSDLC Stage 3: Cross-language test**

    test-geo_pam_cpp_vs_r.R
      expect_equal(geo_pam_cpp(...), geo_pam_r(...), tolerance = 0)
      # Binary matrix — must be EXACT match, not approximate

**Deliverables:** - `src/pam_pip.cpp` — winding number PIP -
`src/pam_rasterize.cpp` — parallel species loop - `src/pam_omp.h` —
thread safety (from xclim) - `src/Makevars`, `src/Makevars.win` —
C++17 + OpenMP - `R/geo_pam.R` — R wrapper with checkmate validation -
`R/geo_pam_r.R` — reference implementation (internal) -
`tests/testthat/test-geo_pam_cpp_vs_r.R`

**Estimated effort:** 2-3 Hermes sessions

### Phase 2: Sparse Storage + S3 Class (v0.4.0)

- `R/pam_class.R` — S3 class with print/summary/plot
- `src/pam_sparse.cpp` — triplet → dgCMatrix conversion
- Richness and range-size extraction methods
- `as_tibble.pam()` for backward compatibility
- Deprecate
  [`get_geo_pam()`](https://alrobles.github.io/pam/reference/get_geo_pam.md)
  with [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html)

**Estimated effort:** 1-2 Hermes sessions

### Phase 3: Environmental PAM (v0.5.0)

**SSDLC Stage 1: R reference**

    env_pam_r(x, env, n_axes, n_bins)
      - terra::extract(env, x) → environmental values
      - prcomp() → PCA scores
      - cut() + interaction() → E-space bins
      - table() → binary matrix

**SSDLC Stage 2: C++ binning engine**

    env_classify_cpp(scores, breaks_per_axis, n_axes)
      - Assign each point to a multidimensional bin
      - Return (bin_id, species_id) triplets

**SSDLC Stage 3: Cross-language test**

**Additional features:** - `auto_resolution()` — optimal bin count via
information-theoretic criterion - Integration with xclim for bioclimatic
variable computation - Support for raw axes (no PCA) or user-supplied
axes

**Estimated effort:** 2-3 Hermes sessions

### Phase 4: Performance + Benchmarks (v0.9.0)

- `inst/benchmarks/bench_geo_pam.R` — benchmark vs terra::rasterize,
  letsR
- `inst/benchmarks/bench_env_pam.R` — benchmark vs ecospat gridding
- Memory profiling (peak RSS for 5000 species × 1° grid)
- Vignette: “Performance Guide” showing speedup curves

**Target benchmarks:** \| Scenario \| v0.2.0 (R) \| v1.0.0 (C++) \|
Speedup \| \|———-\|———–\|————-\|———\| \| 200 spp, 1° \| ~5s \| \<0.5s \|
10× \| \| 5000 spp, 1° \| OOM \| ~30s \| ∞ \| \| 5000 spp, 0.5° \| OOM
\| ~120s \| ∞ \|

### Phase 5: CRAN Prep (v1.0.0)

Following the CRAN Submission Checklist (knowledge base): - \[ \]
`R CMD check --as-cran` → 0 errors, 0 warnings, 0 notes - \[ \] Runnable
`@examples` (not all in `\dontrun`) - \[ \] No `std::cout` — use
`Rcpp::Rcout` only - \[ \] `useDynLib(pam, .registration = TRUE)` in
NAMESPACE - \[ \] `inst/WORDLIST` for proper nouns - \[ \] Vignettes:
“Getting Started”, “Environmental PAM”, “Performance” - \[ \] Test on
win-builder, R-hub, macOS - \[ \] `cran-comments.md`

**Estimated effort:** 1-2 Hermes sessions

------------------------------------------------------------------------

## 7. Delegation Plan for Hermes

Each phase maps to a GitHub issue. Hermes executes on reumanlab/HPC:

    Issue #2: Phase 1 — C++ geographic PAM engine
      Branch: feature/cpp-geo-pam
      Inputs: This design doc + xclim/src/ as reference for Makevars/OpenMP
      Tests: geo_pam_cpp vs geo_pam_r on get_rodentia_mexico() data
      PR target: main

    Issue #3: Phase 2 — Sparse storage + S3 class
      Branch: feature/pam-class
      Depends: Phase 1 merged
      PR target: main

    Issue #4: Phase 3 — Environmental PAM
      Branch: feature/env-pam
      Depends: Phase 2 merged
      PR target: main

    Issue #5: Phase 4 — Benchmarks + optimization
      Branch: feature/benchmarks
      Depends: Phase 3 merged
      PR target: main

    Issue #6: Phase 5 — CRAN submission prep
      Branch: release/v1.0.0
      Depends: Phase 4 merged
      PR target: main

**Hermes task template:**

    Repo: alrobles/PAM
    Branch: feature/<name>
    Design doc: PAM_v1.0.0_DESIGN.md (in repo root)
    Reference code:
      - xclim/src/bioclim_xt.cpp (zero-copy bridge, OpenMP pattern)
      - xclim/src/xclim_omp.h (CRAN-safe thread count)
      - nicher/src/Makevars (RcppParallel/TBB avoidance pattern)
      - maxentcpp/src/rcpp_featured_space.cpp (XPtr + shared_ptr pattern)
    Container: ~/geospatial-rserver/geospatial_latest.sif (R 4.4.2)
    R library: ~/R/rocker-rstudio/4.2
    Test: R CMD check --no-manual --as-cran

------------------------------------------------------------------------

## 8. Scientific Context

### Key References

1.  **Arita, H.T. et al. (2008)** — Species diversity and distribution
    in presence-absence matrices. J. Biogeography. *Theory connecting
    row/column sums to ranges and richness.*

2.  **Arita, H.T. et al. (2012)** — The presence-absence matrix
    reloaded: range-diversity plots. Global Ecol. Biogeography.
    *Numerical procedures for PAM analysis; resolution effects.*

3.  **Soberón, J. & Nakamura, M. (2009)** — Niches and distributional
    areas. PNAS. *Hutchinsonian duality: G-space ↔︎ E-space.*

4.  **Broennimann, O. et al. (2012)** — Measuring ecological niche
    overlap. Ecography. *Environmental-space gridding via PCA.*

5.  **Villalobos, F. et al. (2013)** — Phylogenetic fields of
    species. J. Biogeography. *PAM + phylogeny integration.*

### Relationship to Lab Pipeline

    cochliomyia_hominivorax pipeline:
      01_filter_americas_endemic.R  ← IUCN polygons
      02_rasterize_mammals.R        ← current ad-hoc rasterization
      04_pam_americas_mammals.R     ← reimplements pam logic
                                       (273 lines, manual, slow)

    With pam v1.0.0:
      pam::geo_pam(endemic_vect, res = 1)  ← replaces all of step 04
      pam::env_pam(endemic_vect, bioclim)  ← NEW: environmental PAM

------------------------------------------------------------------------

## 9. Risk Assessment

| Risk | Mitigation |
|----|----|
| CRAN rejects OpenMP | Conditional compilation (xclim pattern), single-thread fallback |
| PIP edge cases (slivers, self-intersecting polygons) | Robust winding-number; pre-validate with [`terra::is.valid()`](https://rspatial.github.io/terra/reference/is.valid.html) |
| Memory for very large PAMs | Sparse storage default; chunked processing for \>10M cells |
| R 4.1.0 minimum too new | PAM targets macroecology researchers on modern R; matches lab standard |
| Boost dependency for R-tree | Use header-only minimal implementation, no external dep |
| Environmental PCA interpretation | Document eigenvalue thresholds; warn if \<70% variance explained |
