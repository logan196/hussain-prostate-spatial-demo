# Live-site audit

Target: https://logan196.github.io/hussain-prostate-spatial-demo/
Method: headless Chrome (Playwright, SwiftShader WebGL) against the **live** GitHub Pages site. Viewports were 1600×950 desktop and 390×844 phone. For each page I recorded console errors, HTTP ≥400 responses, transfer bytes, the time to the first Vitessce canvas, links and horizontal overflow. The script is `ci/smoke.py`.

## Cycle 1: baseline (v0, commit aefe3fd), 2026-10-14 16:10 ET

Pass/fail key: ✗ = problem, ~ = partial, ✓ = OK.

**Severity 1 (fix first)**
1. ✗ **First load: 10.5 MB, 89% of it one file.** `vitessce/he.ome.tif` was a single untiled, uncompressed 9.0 MB strip, downloaded in full before the H&E appeared. `atlas/genes.json` and `atlas/aliases.json` (0.9 MB) also had to load before Vitessce even started.
2. ✗ **No loading state.** The page showed a black panel with one grey line of text for about 3 s, with no sign of progress or of what the page is.
3. ✗ **Phone/narrow widths are broken.** The fixed `1fr 370px` grid gives Vitessce about 20 px at 390 px width. The top bar is clipped and the presets and nav sit off-screen.
4. ✗ **Console errors on every load: 5–6 × 404.** These are `data.zarr/obs/.zarray`, `var/.zarray` and `obs/<column>/.zarray`. The cause is zarr.js probing for an array before falling back to the group. They are harmless, but a scientist who opens DevTools sees red.

**Severity 2**
5. ✗ **State can't be shared.** Presets and answers are not reflected in the URL, so you can't send a colleague "this view".
6. ✗ **No export.** You can't download the answer numbers (CSV) or the current view (PNG).
7. ✗ **The AI features are hard to discover.** The ask box is just a placeholder, and nothing explains LIVE vs PRECOMPUTED or the abstention gating.
8. ✗ **Abstention exists only as a category.** The per-spot voter-agreement feature (`◆ state confidence`) is in the matrix, but no preset shows it as a continuous confidence map.
9. ✗ **No keyboard shortcuts** and no visible focus styles. The preset buttons have no ARIA state and the side panel has no landmark.

**Severity 3**
10. ✗ **The README is stale.** It still describes the old Plotly `index.html` with 337 genes, and it has no GIFs, no citation and no provenance.
11. ✗ **No CI.** Nothing stops a regression from shipping.
12. ~ **Spot Sets panel:** it shows "○" placeholders for the colour swatches under each hierarchy. This is a Vitessce 3.5 rendering quirk; it is cosmetic and left as is.
13. ✓ No broken links (all nav targets returned 200), and no horizontal scroll on desktop. atlas, showcase, brief and explorer had no console errors.

**First impression:** strong. The left Vitessce grid (H&E + UMAP + sets) is good from the first frame, but the first 3 s are blank and the AI layer is invisible until you type something.

## Fix order
| # | Issue | PR |
|---|-------|----|
| 1, 2 | Pyramidal tiled H&E, deferred gene index, boot screen with real progress | #1 |
| 5, 6, 9 | Shareable URLs, CSV/PNG export, keyboard shortcuts, a11y | #3 |
| 3, 7, 8 | Responsive layout, guided tour, confidence-map preset | #4 |
| 10, 11, 4 | README/CHANGELOG, CI smoke test (fails on any console error not on the zarr-probe allowlist) | #5 |

Cycles 2 and 3 are appended below after each merge.
