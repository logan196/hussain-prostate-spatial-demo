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

## Cycle 2: after #1, #3, #4, #5 (live, commit ebd4459), 2026-10-14 ~16:30 ET
- ✓ First-load transfer is **2.9 MB** (was 10.5) on desktop. Vitessce canvases appear at about 3.3 s on a cold load, behind the boot screen.
- ✓ Share URL, CSV and PNG verified on the live site. The PNG is 805×882 with H&E + KLK3 and is not blank. The CSV has 4,394 rows.
- ✓ Phone (390×844): the tissue is visible full-width and there is no horizontal scroll.
- ✗ **First frame still underwhelms**: the tissue filled only about 50% of the spatial view. → fixed by auto-fit in #6.
- ✗ The hypothesis chart axis title overlapped the tick labels. → fixed in #6.
- ✗ The README was stale. → fixed in #7.
- ~ The live smoke run hit one **transient 503** from the GitHub Pages CDN on a zarr chunk, right after a deploy. A re-run passed and `curl` returned 200 three times. This is not a code issue, but CI runs on localhost to avoid it.

## Cycle 3: after #6, #7 (live, commit 1c90237), 2026-10-14 ~16:50 ET
- ✓ `ci/smoke.py` against the **live** site: **9/9 pass**, with zero console errors and zero non-allowlisted HTTP errors.
- ✓ Desktop: 2.9 MB, and the tissue fills the spatial view from the first frame. Atlas, showcase, brief and explorer all have zero errors and no horizontal scroll.

**Remaining issues**
- Phone first load is 7.7 MB, because the fitted zoom pulls full-resolution H&E tiles. A JPEG-compressed pyramid would cut this by about 4×.
- The zarr `.zarray` probe 404s (5 per load) still appear in DevTools. They are inherent to Vitessce's zarr.js reader for AnnData groups and are allowlisted in CI.
- The Vitessce Spot Sets view shows "○" swatch placeholders, and the spatial legend truncates long category names ("Abstained (low co…"). Both are Vitessce rendering behaviours.
- PNG export does not include the Vitessce colour-bar legend, which is an HTML overlay; a caption is burned in instead.
- The GitHub Actions workflow is staged but not active: the CLI token lacks the `workflow` scope.
- The confidence map saturates, because 58% of spots are unanimous (1.0). A discrete 3/5–4/5–5/5 palette would read better.

## Cycle 4: v1.1 pass (live, commits 521b4b4 → 238e047), 2026-10-14 ~17:20 ET
Fixes from the v1.1 plan landed first: JPEG tiles (#9); short labels, discrete confidence levels and the PNG legend (#10). A fresh audit of the live site then found:
1. ✗ **Severity 1: hidden 90 MB download.** The AR vs NE, Gene query and Confidence presets each fired 484 matrix requests. Vitessce's `featureValueHistogram` loads the full obsFeatureMatrix, and this burst also explains the intermittent Pages 503s seen in cycles 2 and 4. → #12 replaces the view with a LIVE Distribution panel and adds a request budget to the smoke test.
2. ✗ **Severity 1 (correctness): proximity answers highlighted the wrong reference set.** Distances are to AR-high spots, but the map showed the ensemble AR call. → #11.
3. ✗ **Severity 2: cell-type questions ("immune cells") returned "unknown".** → #11 adds a labelled marker map.
4. ✗ **Severity 3: the Boundary niche view was off-centre** after the zoom recalibration. → #11.

**Verified live (238e047)**
- `ci/smoke.py`: **11/11 pass**, including the new matrix budget, with zero console errors and zero non-allowlisted HTTP errors.
- First load: 2.1 MB desktop, 2.0 MB phone.

**Remaining**
- The zarr `.zarray` probe 404s are still allowlisted.
- The Actions workflow is still staged only, because the token lacks the `workflow` scope.
- The marker map is one gene per word, so it is a proxy and not a deconvolution.
- The hypothesis chart sits below the fold on 900 px-tall screens.
- The "select all" checkboxes have no tooltip; they use Vitessce's default behaviour.
