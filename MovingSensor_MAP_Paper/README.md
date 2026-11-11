# Moving-sensor inverse source reconstruction — paper experiments

MATLAB code for **“Inverse source problem for the parabolic equation with sparse moving observations”**. This repository provides **only the four source geometries and three experiments per geometry** used in the revised paper: circle, kite, four-leaf, and peanut. No MCMC, alternative timing mode, unrelated demos, or precomputed MAP reconstructions are included.

The forward model uses the authors' polar `Q2` finite elements and backward Euler. Reconstruction uses a constrained Gaussian-prior MAP objective solved by `fmincon` (`sqp`) with one fixed starting point and an optional second SQP refinement. Figures use the paper's original MATLAB drawing routines.

## Requirements

- MATLAB with **Optimization Toolbox** (`fmincon`, `optimoptions`, `sqp`); R2021b is the previously author-used environment.
- Working directory: this repository's top-level folder, containing `CHECK_PROJECT.m`.
- The four `fixtures/V16_REPRO_NOISE_*.mat` files must remain in `fixtures/`; they contain noise **inputs**, never stored reconstructions.
- Sufficient memory for full-order Q2/BE transfer matrices and SVD-based projected optimization.

## Reproduce the paper

```matlab
clear functions
rehash
CHECK_PROJECT

% Run one example from scratch (outputs are NOT loaded from MAT/CSV):
[out,files] = RUN_AND_SAVE_PAPER('PEANUT_EXP2','noisy',1);

% Run all 12 experiment entries from scratch, open paper figures,
% and save new MAT/CSV/TXT files in the ignored results/ directory:
summary = RUN_ALL_AND_SAVE_PAPER('noisy',1);
assert(all([summary.results.completed]));
```

Experiment entry points (call individually as `CIRCLE_EXP1('noisy',1)` etc.):

| Source | Exp. 1: 4 adaptive stages | Exp. 2: 3 fixed sites | Exp. 3: 3 prescribed two-site movements |
|---|---|---|---|
| Circle | `CIRCLE_EXP1` | `CIRCLE_EXP2` | `CIRCLE_EXP3` |
| Kite | `KITE_EXP1` | `KITE_EXP2` | `KITE_EXP3` |
| Four-leaf | `FOURLEAF_EXP1` | `FOURLEAF_EXP2` | `FOURLEAF_EXP3` |
| Peanut | `PEANUT_EXP1` | `PEANUT_EXP2` | `PEANUT_EXP3` |

The outputs are recomputed independently in each experiment. The figures are displayed in MATLAB, **not exported automatically**. Save the relevant Figure windows as EPS if required by the manuscript. Exp. 2 and Exp. 3 create an adaptive reference panel first (`out.figures{1}`); paper comparison panels (a)–(c) correspond to `out.figures{2:4}`.

## One timing rule, with no preassigned future windows

`matlab/V16_ROUTE.m` implements the movement strategy and invokes `V16_NEXT_WINDOW` to choose the **earliest grid-aligned next sensing window** after travel. Only the first sampling time (`1/400`) and the number of samples per window (`80`) are initially set. The current window's noisy endpoint and two noisy auxiliary readings determine the next movement direction and duration. The auxiliary observations are **not** used in MAP fitting. The route is generated causally from measured flux; the MAP fits can be computed afterwards because they do not affect navigation decisions.

With the supplied `('noisy',1)` fixtures, the four observed time windows in units of `1/400` are:

| Geometry | Window 1 | Window 2 | Window 3 | Window 4 |
|---|---|---|---|---|
| Circle | [1,80] | [91,170] | [181,260] | [266,345] |
| Kite | [1,80] | [91,170] | [181,260] | [266,345] |
| Four-leaf | [1,80] | [91,170] | [181,260] | [266,345] |
| **Peanut** | [1,80] | [91,170] | [181,260] | **[265,344]** |

The Peanut fourth window **only** uses grid nodes 265–344; the former 266–345 convention is not an input option in this repository. The expected indices are checked for the fixed manuscript noise fixture, but **not prescribed to the route generator**. `CHECK_PROJECT` also checks the deterministic timing recurrence without computing the PDE.

Each completed four-window experiment has **320 primary observations**. Within each source geometry, the adaptive, fixed-site, and nonadaptive two-site comparisons use identical primary observation times and paired main-measurement Gaussian noise; the same prior/starting point and SQP options are used. The adaptive route additionally consumes six auxiliary flux readings for direction and stopping decisions. Low-flux control sites are deliberately **truth-informed synthetic comparisons**, not a blind site-selection method.

## Results and reproducibility

- `RUN_ALL_AND_SAVE_PAPER` writes a new uniquely named MAT/CSV per experiment and a text log with each example's generated observation windows. Nothing reads old saved MAP results into the solver. `results/` is ignored by Git.
- Each `SAVE_PAPER_RESULT` CSV includes `E_xi`, `E_D`, the MAP objective, and **first-pass / refinement exit flags**. An SQP `exitflag=0` does not mean that a feasible output was unavailable; conversely, feasibility plus objective reduction is not a stationarity certificate.
- The original paper's solver settings, two-pass feasibility-and-improvement acceptance criterion, shape constraints, noise fixtures, figure styling and numerical forward routines are retained. The additional diagnostic fields do not alter the solutions.
- **After changing Peanut's final time grid, its Exp. 1–3 reconstructions, EPS figures and comparison-table entries must be regenerated from this package.** The older reported values from the preassigned-window package are **not** verification targets and are not bundled here.
- This package has undergone archive inspection, cross-file checks and independent timing-arithmetic checks, **but MATLAB is not installed in this packaging environment**, so no claim is made that all 12 new runs or their final parameter/error values have been independently executed here. Please run `CHECK_PROJECT`, the 12 experiments, and inspect EPS figures in MATLAB before publishing figures/tables.

## Publication notes

Before making the repository public, choose the intended software license and replace any placeholder GitHub URL in the article with the actual repository URL. Review any locally generated `results/` text logs for private machine paths before sharing them.
