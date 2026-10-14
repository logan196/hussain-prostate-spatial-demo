# Changelog
The format follows [Keep a Changelog](https://keepachangelog.com/), and the project uses semantic versioning.

## [1.1.0] - 2026-10-14
### Added
- LIVE **Distribution** side-panel section: a log-scale histogram with % detected, median and max for the coloured feature (#12).
- Cell-type words (immune, T cell, macrophage, fibroblast, endothelial, proliferating, …) map to one labelled canonical marker gene (#11).
- **Confidence preset**: 4 discrete voter-agreement levels (5/5, 4/5, 3/5 abstained, 2/5 abstained) plus a Spot Set Sizes chart (#10, #12).
- PNG export now includes a legend: set swatches, or a plasma colour bar over the data min–max (#10).
- Smoke test: a matrix-request budget (fails if a page fires >40 X-chunk requests), plus immune-query and Distribution checks (11 states).

### Changed
- H&E tiles are JPEG q75: 4.8 → 1.27 MB. First load is now 2.1 MB on desktop (was 2.9) and 2.0 MB on phone (was 7.7) (#9).
- **Vitessce's histogram view removed.** It downloaded all 474 matrix chunks (~90 MB, 484 requests) in 3 presets and set off Pages 503s. Every preset now makes 8 matrix requests (#12).
- Short Vitessce category names (`Cell state`: AR adeno, NE / t-SCNC, Amphicrine, Abstained), so legends no longer truncate (#10).
- The Boundary niche preset auto-fits to the niche spots (#11).
- `build_vitessce.py` is synced with the out-of-repo working copy.

### Fixed
- Proximity answers now highlight the AR-high set (AR program z>1, n=628) that distances are measured to, instead of the ensemble AR set (n=515) (#11).
- The Spot Sets "select all" checkboxes no longer wrap onto their own line as bare "○" (#10).

## [1.0.0] - 2026-10-14
Three audit → fix → verify cycles against the live Pages site. See `docs/AUDIT.md`.

### Added
- Branded boot screen whose progress comes from real resource timings (#1).
- Shareable view URLs, per-spot CSV export with a provenance header, spatial PNG export, keyboard shortcuts and a help dialog (#3).
- **Confidence** preset, a per-spot voter-agreement map of the abstention-gated ensemble, plus a 4-step guided tour (#4).
- Responsive layout for phones and narrow windows (#4).
- `ci/lint.mjs` and a `ci/smoke.py` headless smoke test that fails on console errors, with the Actions workflow staged at `ci/github-workflow-ci.yml` (#5).
- `docs/AUDIT.md` (#2), README with GIFs, citation and provenance, and this changelog (#7).

### Changed
- H&E is now a pyramidal, tiled, zlib OME-TIFF (9.0 → 4.8 MB on disk, and only visible tiles are fetched). First-load transfer fell from 10.5 MB to 2.9 MB (#1).
- The gene index and aliases load after first paint (#1).
- Presets auto-fit the tissue to the spatial view (#6).

### Fixed
- Hypothesis chart axis title overlapping the tick labels (#6).
- "compare <gene>" queries now colour the views by that gene.

## [0.9.0] - 2026-10-14
- First public release: Vitessce core with NL query, abstention-gated ensemble, insight and hypothesis panels; gene atlas; deck.gl showcase; brief and plan PDF.
