# Hussain prostate spatial demo: AR → NE lineage in space

An interactive spatial explorer and a reproducible analysis of **GEO GSE230282**: Visium CytAssist FFPE prostate with co-existing de novo NEPC and ARPC (Watanabe et al., *IJMS* 2023). Prepared for a call with Dr. Maha Hussain (Northwestern) about treatment-emergent small-cell NE prostate cancer.

- `index.html` + `data/`: static portal (Plotly; runs on GitHub Pages, no server). Pan/zoom the H&E, toggle program, score and cluster overlays, pick any of 337 genes, and box/lasso-select spots for summary stats.
- `analysis.py`: re-runs everything; downloads the data from GEO if it is missing. `pip install scanpy squidpy scikit-misc leidenalg igraph`, then `python analysis.py [--z 1.0]`.
- `stats.json`: every number quoted in the brief.
- `figures/`: static PNGs.
- `brief.md` / `brief.html` and `hussain_profile.md`: call prep.

Portal styling is reused from Logan's `advocate-health-portal` (Novaflow IDE theme).

## Gene atlas (atlas.html)
Gene-centric atlas of all 18,939 QC-passing genes. Hash routes: `#/` is the searchable index and `#/gene/<SYMBOL>` is the gene page, e.g. `#/gene/KLK3`. Optional query parameters: `?mode=region|cl|he` and `?zoom=N`.
Each gene page has:
- HGNC metadata, with links to NCBI, GeneCards, HPA and Ensembl.
- A zoomable H&E with an expression/region/cluster overlay.
- Expression by region and by cluster.
- NE vs AR Wilcoxon statistics.
- The top 20 correlated genes.
- Provenance.

Build the atlas with `python build_atlas.py` (after `analysis.py`). Data lives in `atlas/`: uint8 expression shards, per-gene detail JSON, `genes.json` and `meta.json`. The earlier single-page viewer is now `explorer.html`.

## 3D & animations (showcase.html)
- **showcase.html** (deck.gl 9, WebGL2). Four scenes:
  - A dual-channel AR (blue) / NE (magenta) glow map, blended additively so co-expression shows as violet.
  - A 3D extruded AR/NE landscape that grows in on load.
  - An animated AR→NE sweep with a live NE/AR-vs-distance chart.
  - A spatial↔UMAP morph.

  It also has camera fly-to presets, a ▶ Tour fly-through, and animated stat counters and insight cards. URL parameters: `?scene=glow|land|sweep|umap&view=overview|tilt|ne|ar|edge&tour=1`.
- **vitessce.html** (Vitessce 3.5.1 from CDN). Linked views: spatial (H&E OME-TIFF + spots), UMAP, gene list, region/Leiden sets, layer controller and expression histogram. Data is `vitessce/data.zarr` (327 genes, AnnData zarr v2) plus `vitessce/he.ome.tif`. It needs HTTP Range support, which GitHub Pages provides.

## Core app: Vitessce (index.html)
A full-screen Vitessce 3.5.1 app. Slim top bar with four presets: Overview / AR vs NE / Boundary niche / Gene query.
- **Data:** `vitessce/data.zarr` holds all 18,939 genes plus 5 derived ◆ features, as AnnData zarr v2 in 40-gene chunks with blosc-zstd compression (136 MB).
- **H&E:** `vitessce/he.ome.tif` (planar CYX) is displayed in true color via `photometricInterpretation: {init_A_image_0: "RGB"}`.
- **Ask box (LIVE):** parsed client-side into intent, gene (exact, HGNC alias, or Levenshtein fuzzy match), region and distance. The parse is applied as Vitessce coordination changes: feature selection, color encoding, and new spot sets (additionalObsSets). Answers are computed in the browser from `atlas/expr` shards and `vitessce/spots.json`.
- **Ensemble cell-state calls (PRECOMPUTED, build_vitessce.py):** 5 voters (z>1 thresholds, relaxed z>0.75, spatially smoothed, 8-component GMM, Leiden-cluster profile) vote on AR / NE / Amphicrine / Stroma / Other. A spot abstains unless at least 4 of 5 voters agree.
- **Insight cards (PRECOMPUTED):** generated per preset from `stats.json` and the ensemble output.
- **Hussain hypothesis panel (PRECOMPUTED):** distance gradient, Spearman correlation, and transition-niche numbers.
- Deep links: `?preset=overview|arne|niche|gene` and `?q=<question>`.
