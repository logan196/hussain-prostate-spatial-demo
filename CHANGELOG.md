# Changelog
The format follows [Keep a Changelog](https://keepachangelog.com/), and the project uses semantic versioning.

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
