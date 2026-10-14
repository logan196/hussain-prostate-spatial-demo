# Hussain prostate spatial demo: AR → NE lineage in space

An interactive spatial explorer and a reproducible analysis of **GEO GSE230282**: Visium CytAssist FFPE prostate with co-existing de novo NEPC and ARPC (Watanabe et al., *IJMS* 2023). Prepared for a call with Dr. Maha Hussain (Northwestern) about treatment-emergent small-cell NE prostate cancer.

- `index.html` + `data/`: static portal (Plotly; runs on GitHub Pages, no server). Pan/zoom the H&E, toggle program, score and cluster overlays, pick any of 337 genes, and box/lasso-select spots for summary stats.
- `analysis.py`: re-runs everything; downloads the data from GEO if it is missing. `pip install scanpy squidpy scikit-misc leidenalg igraph`, then `python analysis.py [--z 1.0]`.
- `stats.json`: every number quoted in the brief.
- `figures/`: static PNGs.
- `brief.md` / `brief.html` and `hussain_profile.md`: call prep.

Portal styling is reused from Logan's `advocate-health-portal` (Novaflow IDE theme).

## Gene atlas (index.html)
Gene-centric atlas of all 18,939 QC-passing genes. Hash routes: `#/` is the searchable index and `#/gene/<SYMBOL>` is the gene page, e.g. `#/gene/KLK3`. Optional query parameters: `?mode=region|cl|he` and `?zoom=N`.
Each gene page has:
- HGNC metadata, with links to NCBI, GeneCards, HPA and Ensembl.
- A zoomable H&E with an expression/region/cluster overlay.
- Expression by region and by cluster.
- NE vs AR Wilcoxon statistics.
- The top 20 correlated genes.
- Provenance.

Build the atlas with `python build_atlas.py` (after `analysis.py`). Data lives in `atlas/`: uint8 expression shards, per-gene detail JSON, `genes.json` and `meta.json`. The earlier single-page viewer is now `explorer.html`.
