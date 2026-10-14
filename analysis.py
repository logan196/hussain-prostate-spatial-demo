#!/usr/bin/env python3
"""AR vs neuroendocrine (NE) spatial analysis of GSE230282 (Visium CytAssist, FFPE prostate,
co-existing ARPC + NEPC; Watanabe et al., IJMS 2023, PMID 37240308).
Usage: python analysis.py [--data data] [--out .]   (downloads GEO files if missing)"""
import argparse, gzip, json, os, shutil, urllib.request
import numpy as np, pandas as pd, scanpy as sc, squidpy as sq, scipy.io, scipy.sparse as sp
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7211nnn/GSM7211257/suppl/GSM7211257_EHU-W3_"
FILES = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz", "scalefactors_json.json.gz",
         "tissue_hires_image.png.gz", "tissue_positions.csv.gz"]
AR = ["AR", "KLK3", "KLK2", "TMPRSS2", "NKX3-1", "FKBP5", "STEAP2", "PMEPA1", "SLC45A3", "HOXB13", "FOLH1"]
NE = ["CHGA", "CHGB", "SYP", "ASCL1", "INSM1", "NCAM1", "ENO2", "SOX2", "SCG2", "SCG3", "PCSK1", "STMN2", "ELAVL4", "NEUROD1"]
EXTRA = ["MKI67", "TOP2A", "RB1", "TP53", "PTEN", "EZH2", "MYCN", "AURKA", "FOXA2", "POU2F3", "YAP1", "DLL3",
         "SOX9", "KRT8", "KRT18", "KRT5", "TP63", "VIM", "ACTA2", "COL1A1", "PTPRC", "CD68", "PECAM1",
         "HGF", "HMOX1", "ELN", "GREM1", "RBFOX3", "SFTPC"]
UM_PER_SPOT_DIAM = 55.0

def fetch(d):
    os.makedirs(d, exist_ok=True)
    for f in FILES:
        p = os.path.join(d, "GSM7211257_EHU-W3_" + f)
        if not os.path.exists(p):
            print("downloading", f); urllib.request.urlretrieve(GEO + f, p)
    return lambda f: os.path.join(d, "GSM7211257_EHU-W3_" + f)

def load(P):
    X = scipy.io.mmread(gzip.open(P("matrix.mtx.gz"))).T.tocsr()
    bc = pd.read_csv(P("barcodes.tsv.gz"), header=None)[0].values
    ft = pd.read_csv(P("features.tsv.gz"), header=None, sep="\t")
    ad = sc.AnnData(X.astype(np.float32), obs=pd.DataFrame(index=bc), var=pd.DataFrame({"gene_id": ft[0].values}, index=ft[1].values))
    ad.var_names_make_unique()
    pos = pd.read_csv(P("tissue_positions.csv.gz"), index_col=0)
    ad = ad[ad.obs_names.isin(pos.index[pos.in_tissue == 1])].copy()
    pos = pos.loc[ad.obs_names]
    ad.obs[["array_row", "array_col"]] = pos[["array_row", "array_col"]].values
    ad.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    sf = json.load(gzip.open(P("scalefactors_json.json.gz")))
    img = np.array(Image.open(gzip.open(P("tissue_hires_image.png.gz"))).convert("RGB"))
    ad.uns["spatial"] = {"EHU-W3": {"images": {"hires": img}, "scalefactors": sf}}
    return ad, sf, img

def main():
    a = argparse.ArgumentParser(); a.add_argument("--data", default="data"); a.add_argument("--out", default=".")
    a.add_argument("--z", type=float, default=1.0, help="z-score threshold for calling a program 'high'")
    args = a.parse_args(); O = lambda *p: os.path.join(args.out, *p)
    os.makedirs(O("figures"), exist_ok=True); os.makedirs(O("portal", "data"), exist_ok=True)
    ad, sf, img = load(fetch(args.data))
    um_px = UM_PER_SPOT_DIAM / sf["spot_diameter_fullres"]
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    n0 = ad.n_obs; sc.pp.filter_cells(ad, min_counts=500); sc.pp.filter_genes(ad, min_cells=10)
    ad.layers["counts"] = ad.X.copy(); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    arg = [g for g in AR if g in ad.var_names]; neg = [g for g in NE if g in ad.var_names]
    sc.tl.score_genes(ad, arg, score_name="AR_score", random_state=0)
    sc.tl.score_genes(ad, neg, score_name="NE_score", random_state=0)
    z = lambda v: (v - v.mean()) / v.std()
    ad.obs["AR_z"] = z(ad.obs.AR_score); ad.obs["NE_z"] = z(ad.obs.NE_score)
    arh = ad.obs.AR_z > args.z; neh = ad.obs.NE_z > args.z
    cls = np.where(arh & neh, "Dual AR+/NE+", np.where(arh, "AR-high", np.where(neh, "NE-high", "Low/other")))
    ad.obs["program"] = pd.Categorical(cls, categories=["AR-high", "NE-high", "Dual AR+/NE+", "Low/other"])
    ad.obs["NE_minus_AR"] = ad.obs.NE_z - ad.obs.AR_z
    # clustering
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.pca(ad, n_comps=30, mask_var="highly_variable", random_state=0); sc.pp.neighbors(ad, random_state=0)
    sc.tl.leiden(ad, resolution=0.5, random_state=0, flavor="igraph", n_iterations=2, directed=False); sc.tl.umap(ad, random_state=0)
    # distances (um)
    xy = ad.obsm["spatial"] * um_px
    ar_xy, ne_xy = xy[arh.values & ~neh.values], xy[neh.values & ~arh.values]
    d_ne, _ = cKDTree(ar_xy).query(ne_xy); d_rand = cKDTree(ar_xy).query(xy[~arh.values])[0]
    ad.obs["dist_to_AR_um"] = cKDTree(ar_xy).query(xy)[0]
    # neighborhood
    sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
    sq.gr.nhood_enrichment(ad, cluster_key="program", seed=0, n_perms=1000)
    A = ad.obsp["spatial_connectivities"]
    nb_ne = np.asarray(A @ neh.values.astype(float)).ravel() / np.maximum(np.asarray(A.sum(1)).ravel(), 1)
    ad.obs["frac_NE_neighbors"] = nb_ne
    # border (AR spots touching NE spots) vs interior
    pure_ar = (arh & ~neh).values; pure_ne = (neh & ~arh).values
    border_ar = pure_ar & (np.asarray(A @ pure_ne.astype(float)).ravel() > 0)
    corr = np.corrcoef(ad.obs.AR_score, ad.obs.NE_score)[0, 1]
    vc = ad.obs.program.value_counts()
    bins = [0, 100, 250, 500, 1000, 2000, 1e9]; lab = ["0-100", "100-250", "250-500", "500-1000", "1000-2000", ">2000"]
    gradient = ad.obs.assign(b=pd.cut(ad.obs.dist_to_AR_um, bins, labels=lab, include_lowest=True)).groupby("b", observed=False)[["NE_z", "AR_z"]].mean().round(3)
    cl_tab = ad.obs.groupby("leiden", observed=True)[["AR_z", "NE_z"]].mean().round(2).assign(n=ad.obs.leiden.value_counts().sort_index())
    Z = ad.uns["program_nhood_enrichment"]["zscore"]
    stats = {
        "dataset": "GSE230282 / GSM7211257 (EHU-W3), Visium CytAssist FFPE",
        "spots_in_tissue": int(n0), "spots_after_qc": int(ad.n_obs), "genes_after_qc": int(ad.n_vars),
        "median_umi": float(np.median(ad.obs.total_counts)), "median_genes": float(np.median(ad.obs.n_genes_by_counts)),
        "AR_genes_used": arg, "NE_genes_used": neg, "z_threshold": args.z,
        "program_counts": {k: int(v) for k, v in vc.items()},
        "program_fraction_pct": {k: round(100 * v / ad.n_obs, 2) for k, v in vc.items()},
        "dual_fraction_of_all_spots_pct": round(100 * (arh & neh).mean(), 2),
        "dual_fraction_of_program_high_spots_pct": round(100 * (arh & neh).sum() / (arh | neh).sum(), 2),
        "AR_NE_score_pearson_r": round(float(corr), 3),
        "NE_to_nearest_AR_um": {"median": round(float(np.median(d_ne)), 1), "mean": round(float(d_ne.mean()), 1),
                                "p25": round(float(np.percentile(d_ne, 25)), 1), "p75": round(float(np.percentile(d_ne, 75)), 1),
                                "max": round(float(d_ne.max()), 1), "pct_within_250um": round(100 * float((d_ne <= 250).mean()), 1)},
        "nonAR_to_nearest_AR_um_median": round(float(np.median(d_rand)), 1),
        "AR_spots_bordering_NE": int(border_ar.sum()), "AR_spots_total_pure": int(pure_ar.sum()),
        "nhood_enrichment_z": {r: {c: round(float(Z[i, j]), 1) for j, c in enumerate(ad.obs.program.cat.categories)} for i, r in enumerate(ad.obs.program.cat.categories)},
        "NE_AR_by_distance_to_AR_um": gradient.reset_index().rename(columns={"b": "bin"}).to_dict("records"),
        "leiden_clusters": cl_tab.reset_index().to_dict("records"),
        "um_per_fullres_px": um_px,
    }
    json.dump(stats, open(O("stats.json"), "w"), indent=1, default=str); print(json.dumps(stats, indent=1, default=str))
    # figures
    kw = dict(library_id="EHU-W3", img_key="hires", size=1.4, show=False, frameon=False)
    for k, cm, t in [("AR_score", "Reds", "AR program score"), ("NE_score", "Blues", "NE program score")]:
        sc.pl.spatial(ad, color=k, cmap=cm, title=t, **kw); plt.savefig(O("figures", f"{k}.png"), dpi=160, bbox_inches="tight"); plt.close()
    pal = {"AR-high": "#e4572e", "NE-high": "#2e86de", "Dual AR+/NE+": "#a855f7", "Low/other": "#d0d0d0"}
    ad.uns["program_colors"] = [pal[c] for c in ad.obs.program.cat.categories]
    sc.pl.spatial(ad, color="program", title="Program calls (z>%.1f)" % args.z, **kw); plt.savefig(O("figures", "program_calls.png"), dpi=160, bbox_inches="tight"); plt.close()
    sc.pl.spatial(ad, color="leiden", title="Leiden clusters", **kw); plt.savefig(O("figures", "leiden.png"), dpi=160, bbox_inches="tight"); plt.close()
    f, ax = plt.subplots(figsize=(5, 4.5)); ax.scatter(ad.obs.AR_z, ad.obs.NE_z, s=3, c=[pal[c] for c in ad.obs.program], alpha=.6)
    ax.axvline(args.z, ls="--", c="k", lw=.6); ax.axhline(args.z, ls="--", c="k", lw=.6); ax.set_xlabel("AR score (z)"); ax.set_ylabel("NE score (z)")
    ax.set_title(f"AR vs NE per spot (r={corr:.2f})"); f.savefig(O("figures", "AR_vs_NE_scatter.png"), dpi=160, bbox_inches="tight"); plt.close()
    f, ax = plt.subplots(figsize=(5.5, 3.5)); ax.hist(d_ne, bins=40, color="#2e86de"); ax.set_xlabel("NE-high spot → nearest AR-high spot (µm)"); ax.set_ylabel("spots")
    ax.set_title(f"median {np.median(d_ne):.0f} µm"); f.savefig(O("figures", "NE_to_AR_distance.png"), dpi=160, bbox_inches="tight"); plt.close()
    sq.pl.nhood_enrichment(ad, cluster_key="program", title="Neighborhood enrichment (z)", cmap="coolwarm", vmin=-50, vmax=50); plt.savefig(O("figures", "nhood_enrichment.png"), dpi=160, bbox_inches="tight"); plt.close()
    f, ax = plt.subplots(figsize=(5.5, 3.5)); gradient.plot(ax=ax, marker="o", color=["#2e86de", "#e4572e"]); ax.set_xlabel("distance to nearest AR-high spot (µm)"); ax.set_ylabel("mean z")
    ax.set_title("NE/AR gradient away from AR regions"); f.savefig(O("figures", "gradient.png"), dpi=160, bbox_inches="tight"); plt.close()
    sc.pl.dotplot(ad, arg + neg, groupby="leiden", show=False, standard_scale="var"); plt.savefig(O("figures", "dotplot_leiden.png"), dpi=160, bbox_inches="tight"); plt.close()
    # portal export
    s = sf["tissue_hires_scalef"]; im = Image.fromarray(img); W = 1400; r = W / im.width
    im.resize((W, int(im.height * r))).save(O("portal", "data", "tissue.jpg"), quality=85)
    xyh = ad.obsm["spatial"] * s * r
    genes = [g for g in dict.fromkeys(arg + neg + EXTRA) if g in ad.var_names]
    hv = ad.var[ad.var.highly_variable].sort_values("highly_variable_rank").index[:300].tolist()
    genes = list(dict.fromkeys(genes + hv))
    M = ad[:, genes].X; M = M.toarray() if sp.issparse(M) else M
    out = {"w": W, "h": int(im.height * r), "spot_px": sf["spot_diameter_fullres"] * s * r, "um_per_px": um_px / (s * r),
           "x": np.round(xyh[:, 0], 1).tolist(), "y": np.round(xyh[:, 1], 1).tolist(), "barcode": ad.obs_names.tolist(),
           "obs": {k: np.round(ad.obs[k].values.astype(float), 3).tolist() for k in ["AR_score", "NE_score", "AR_z", "NE_z", "NE_minus_AR", "dist_to_AR_um", "frac_NE_neighbors", "total_counts"]},
           "program": ad.obs.program.astype(str).tolist(), "leiden": ad.obs.leiden.astype(str).tolist(),
           "umap": np.round(ad.obsm["X_umap"], 2).tolist(), "genes": genes, "gene_sets": {"AR": arg, "NE": neg},
           "expr": {g: np.round(M[:, i], 2).tolist() for i, g in enumerate(genes)}, "stats": stats}
    json.dump(out, open(O("portal", "data", "spots.json"), "w"), separators=(",", ":"), default=str)
    ad.write_h5ad(O("data", "processed.h5ad"))

if __name__ == "__main__":
    main()
