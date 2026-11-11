# Release changes and validation scope

This release intentionally contains **one** implementation of the 12 manuscript experiments, rather than the original archive's separate fixed-window paper drivers and online-timing demo.

## Changes relative to the supplied 2026-11-09 ZIP

1. `matlab/V16_ROUTE.m` now uses the original package's tested recursive online scheduling logic to construct every successive observation window from a measured-flux decision and main-sensor travel time. `matlab/V16_NEXT_WINDOW.m` is the shared grid rounding helper. The 12 `*_EXP*.m` entry points, meshes, MAP objective and figure styles are unchanged.
2. `matlab/V16_CONFIG.m` removes the preassigned future `cfg.windows`; the only initial time is `cfg.initialStart=cfg.dt`. `MAP_PREPARE_CASE` verifies the fixed supplied-noise run has four windows and the expected positions/times. The expected time values are **assertions after route generation**, not values used in routing.
3. `matlab/V16_BUILD_DATA.m` constructs all seven comparison policies from the single generated time grid. Peanut's fourth window is `[265,344]/400`, for adaptive, all fixed controls and all two-site controls.
4. `matlab/solve_map_sqp.m` records the refinement SQP exit flag and whether the refinement was accepted. `SAVE_PAPER_RESULT.m` exports these diagnostics. **The optimizer settings and acceptance test were not changed**.
5. Removed old timing demos/audits, extra scripts and documents, and the stale `reference_results/` CSVs produced from the earlier Peanut time grid. Neither the old preassigned-window timing path nor its archived outputs remain runnable from this release.

## What was and was not validated here

- Static inspection verified the sole production route and policy builder have no dependency on `cfg.windows` or the removed online-timing entry points.
- Deterministic timing arithmetic was independently checked for all shapes: stage indices 1:80, 91:170, 181:260, and fourth stage 266:345 (circle/kite/four-leaf) or 265:344 (peanut).
- Noise fixtures and original Q2/BE routines are kept unchanged. Output file integrity and direct-source dependencies were reviewed.
- **No MATLAB or Optimization Toolbox is available in this packaging container.** All 12 revised output estimates have **not** been re-run or numerically verified here. `README.md` describes the exact MATLAB commands required to regenerate the manuscript figures and tables.
- The supplied prior author log confirms the independent dynamic Peanut route's fourth window [265,344]/400 and 320 primary readings. This evidence is not a substitute for re-running the three main-paper Peanut comparisons.
