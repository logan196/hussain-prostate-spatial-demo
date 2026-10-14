# Hussain prostate spatial demo: AR → NE lineage in space

An interactive spatial explorer and a reproducible analysis of **GEO GSE230282**: Visium CytAssist FFPE prostate with co-existing de novo NEPC and ARPC (Watanabe et al., *IJMS* 2023). Prepared for a call with Dr. Maha Hussain (Northwestern) about treatment-emergent small-cell NE prostate cancer.

- `index.html` + `data/`: static portal (Plotly; runs on GitHub Pages, no server). Pan/zoom the H&E, toggle program, score and cluster overlays, pick any of 337 genes, and box/lasso-select spots for summary stats.
- `analysis.py`: re-runs everything; downloads the data from GEO if it is missing. `pip install scanpy squidpy scikit-misc leidenalg igraph`, then `python analysis.py [--z 1.0]`.
- `stats.json`: every number quoted in the brief.
- `figures/`: static PNGs.
- `brief.md` / `brief.html` and `hussain_profile.md`: call prep.

Portal styling is reused from Logan's `advocate-health-portal` (Novaflow IDE theme).
