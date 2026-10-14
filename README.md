# Prostate Spatial Atlas: AR → NE lineage in space

**Live:** https://logan196.github.io/hussain-prostate-spatial-demo/

This is an interactive spatial-transcriptomics portal for **GEO GSE230282 / GSM7211257**: one Visium CytAssist FFPE prostate sample in which de novo neuroendocrine prostate cancer (NEPC) sits alongside AR-driven adenocarcinoma (Watanabe et al., *IJMS* 2023). It was built as call prep for Dr. Maha Hussain (Northwestern) on treatment-emergent small-cell NE prostate cancer (t-SCNC).

![Vitessce core: presets, confidence map and a live NL query](docs/media/vitessce_tour.gif)

## What's in it
| Page | What it is |
|---|---|
| `index.html` | **Vitessce 3.5 core.** H&E (pyramidal OME-TIFF) + UMAP + spot sets + all 18,939 genes, lazily loaded per 40-gene zarr chunk. It has 5 presets (Overview, AR vs NE, Boundary niche, Gene query, Confidence) and a natural-language ask box. The right panel holds LIVE answers and PRECOMPUTED insights, the hypothesis panel and the ensemble cell-state calls. |
| `atlas.html` | Gene atlas, one page per gene (`#/gene/KLK3`): HGNC metadata, a zoomable overlay, expression by region and cluster, NE-vs-AR Wilcoxon stats, top correlates and provenance. |
| `showcase.html` | deck.gl 3D and animations: AR/NE glow map, extruded landscape, AR→NE sweep, spatial↔UMAP morph. |
| `brief.html`, `Novaflow_Hussain_Plan.pdf` | Call brief and the 4-page plan with references. |
| `explorer.html` | The earlier Plotly viewer (337 genes), kept for comparison. |

![3D landscape (showcase.html)](docs/media/anim_landscape_grow.gif)

### Using the core viewer
- **Ask**: e.g. `compare KLK3 in AR vs NE`, `NE spots within 300 µm of AR`, `where is chromogranin high` (aliases and fuzzy symbols work), or `genes correlated with SYP`. The parser runs in your browser and is rule-based: no LLM and no server. Each answer is computed from per-spot data and shown with the intent it parsed.
- **Share**: the URL always encodes the view (`?preset=conf`, `?q=compare%20KLK3%20in%20AR%20vs%20NE`, `?preset=gene&gene=SYP`). The 🔗 button copies it.
- **Export**:
  - CSV of every spot: ensemble call, voter agreement, distance to AR, region, the current gene and the query group, with a provenance header.
  - PNG of the spatial view.
- **Keys**: `1`–`5` switch presets, `/` focuses the ask box, `S` copies the link, `E` exports CSV, `P` exports PNG, `?` shows help. The tour is under **Tour** (or `?tour=1`).

## Honesty rules
- **LIVE** = computed in the browser for your query. **PRECOMPUTED** = taken from `stats.json` and the ensemble produced by `analysis.py`.
- **Abstention gating**: 5 voters (thresholds, relaxed thresholds, spatial smoothing, GMM, cluster profile) each call every spot. A label is kept only if ≥4/5 agree; otherwise the spot is **Abstained**. That happens for 889 spots (20.23%). The **Confidence** preset maps per-spot agreement.
- **Limits**: this is one patient with de novo NEPC and no treatment history. The treatment-sequence → t-SCNC hypothesis is shown as **not testable** with these data. All statistics are spot-level and descriptive.

## Reproduce
```bash
pip install scanpy squidpy scikit-misc leidenalg igraph anndata zarr tifffile
python analysis.py          # downloads GSE230282 from GEO if missing → stats.json, figures/
python build_atlas.py       # atlas/ shards
python build_vitessce.py    # vitessce/data.zarr, pyramidal he.ome.tif, spots.json
python -m RangeHTTPServer 8765   # the OME-TIFF needs HTTP Range support
```

## Quality checks
- `node ci/lint.mjs`: syntax-checks every inline script and runs htmlhint.
- `python ci/smoke.py http://localhost:8765/`: drives headless Chromium through 9 page states and **fails on any console error or HTTP ≥400**. The only allowlisted errors are zarr.js `.zarray` group probes.
- The GitHub Actions workflow is in `ci/github-workflow-ci.yml`. To enable it, move it to `.github/workflows/ci.yml` (pushing workflow files requires a token with the `workflow` scope).
- The audit log is in [`docs/AUDIT.md`](docs/AUDIT.md), and changes are in [`CHANGELOG.md`](CHANGELOG.md).

## Provenance & citation
- **Data**: GEO [GSE230282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230282), sample GSM7211257 (10x Visium CytAssist, FFPE).
- **Primary paper**: Watanabe R, et al. *Spatial gene expression analysis reveals characteristic gene expression patterns of de novo neuroendocrine prostate cancer coexisting with androgen receptor pathway prostate cancer.* **Int J Mol Sci** 2023;24(10):8955. doi:[10.3390/ijms24108955](https://doi.org/10.3390/ijms24108955), PMID [37240308](https://pubmed.ncbi.nlm.nih.gov/37240308/).
- **Software**:
  - Vitessce (Keller et al., *Nat Methods* 2025) and deck.gl.
  - scanpy (Wolf et al., *Genome Biol* 2018) and squidpy (Palla et al., *Nat Methods* 2022).
  - Gene metadata from HGNC.
- **Derived values**: everything shown is computed from the public data by the scripts in this repo. Nothing is simulated.

Please cite the primary paper and GEO accession if you reuse the data. Cite this repo as: Beharry L. *Prostate Spatial Atlas (GSE230282 demo)*, v1.0, 2026. https://github.com/logan196/hussain-prostate-spatial-demo

Portal styling is adapted from Logan's `advocate-health-portal` (Novaflow theme).
