#!/usr/bin/env python3
"""Build static gene-centric atlas data (portal/atlas/) from data/processed.h5ad (output of analysis.py)."""
import json, os, subprocess, hashlib, gzip
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy.spatial import cKDTree
from PIL import Image
OUT = "portal/atlas"; SHARD = 200; TOPC = 20
os.makedirs(OUT + "/expr", exist_ok=True); os.makedirs(OUT + "/detail", exist_ok=True)
ad = sc.read_h5ad("data/processed.h5ad")
sf = json.load(gzip.open("data/GSM7211257_EHU-W3_scalefactors_json.json.gz"))
img = Image.open(gzip.open("data/GSM7211257_EHU-W3_tissue_hires_image.png.gz")).convert("RGB")
img.save(OUT + "/tissue_hires.jpg", quality=92)
s = sf["tissue_hires_scalef"]; xy = ad.obsm["spatial"] * s
um_px = 55.0 / sf["spot_diameter_fullres"]
# regions
prog = ad.obs.program.astype(str).values
A = ad.obsp["spatial_connectivities"]
isA = prog == "AR-high"; isN = prog == "NE-high"
nbA = np.asarray(A @ isA.astype(float)).ravel() > 0; nbN = np.asarray(A @ isN.astype(float)).ravel() > 0
low = prog == "Low/other"
region = np.where(isA, "AR region", np.where(isN, "NE region", np.where(prog == "Dual AR+/NE+", "Dual AR+/NE+",
          np.where(low & nbA, "AR boundary", np.where(low & nbN, "NE boundary", "Stroma/other")))))
REG = ["AR region", "AR boundary", "Dual AR+/NE+", "Stroma/other", "NE boundary", "NE region"]
ad.obs["region"] = pd.Categorical(region, categories=REG)
cl = ad.obs.leiden.astype(str).values; CL = sorted(set(cl), key=int)
X = ad.X.toarray().astype(np.float32) if sp.issparse(ad.X) else np.asarray(ad.X, np.float32)
genes = ad.var_names.tolist(); G = len(genes)
# DE NE vs AR
sc.tl.rank_genes_groups(ad, "program", groups=["NE-high"], reference="AR-high", method="wilcoxon", pts=False)
r = ad.uns["rank_genes_groups"]
de = pd.DataFrame({"g": r["names"]["NE-high"], "lfc": r["logfoldchanges"]["NE-high"], "padj": r["pvals_adj"]["NE-high"]}).set_index("g")
# region/cluster stats
def grp(lab, cats):
    m = np.zeros((len(cats), G), np.float32); p = np.zeros_like(m); n = []
    for i, c in enumerate(cats):
        k = lab == c; n.append(int(k.sum())); m[i] = X[k].mean(0); p[i] = (X[k] > 0).mean(0)
    return m, p, n
rm, rp, rn = grp(region, REG); cm, cp, cn = grp(cl, CL)
mean_all = X.mean(0); pct_all = (X > 0).mean(0); mx = X.max(0)
# correlations (blockwise, z-scored)
Z = (X - mean_all) / np.maximum(X.std(0), 1e-6); Z /= np.sqrt(X.shape[0])
cor_ar = Z.T @ ((ad.obs.AR_score.values - ad.obs.AR_score.mean()) / ad.obs.AR_score.std() / np.sqrt(X.shape[0])).astype(np.float32)
cor_ne = Z.T @ ((ad.obs.NE_score.values - ad.obs.NE_score.mean()) / ad.obs.NE_score.std() / np.sqrt(X.shape[0])).astype(np.float32)
topi = np.zeros((G, TOPC), np.int32); topv = np.zeros((G, TOPC), np.float32)
for b in range(0, G, 500):
    C = Z[:, b:b+500].T @ Z
    for j in range(C.shape[0]): C[j, b+j] = -9
    C[:, pct_all < 0.05] = -9  # partners must be detected in >=5% of spots
    idx = np.argpartition(-C, TOPC, axis=1)[:, :TOPC]
    v = np.take_along_axis(C, idx, 1); o = np.argsort(-v, 1)
    topi[b:b+500] = np.take_along_axis(idx, o, 1); topv[b:b+500] = np.take_along_axis(v, o, 1)
del Z
# HGNC metadata
h = pd.read_csv("data/hgnc.tsv", sep="\t", usecols=["symbol", "name", "entrez_id", "ensembl_gene_id", "location", "locus_group"], dtype=str)
hs = h.set_index("symbol"); he = h.dropna(subset=["ensembl_gene_id"]).set_index("ensembl_gene_id")
AR = set(ad.uns.get("AR_genes", [])) or {"AR","KLK3","KLK2","TMPRSS2","NKX3-1","FKBP5","STEAP2","PMEPA1","SLC45A3","HOXB13","FOLH1"}
NE = {"CHGA","CHGB","SYP","ASCL1","INSM1","NCAM1","ENO2","SOX2","SCG2","SCG3","PCSK1","STMN2","ELAVL4","NEUROD1"}
index = []
for gi, g in enumerate(genes):
    ens = ad.var.gene_id.iloc[gi]
    row = he.loc[ens] if ens in he.index else (hs.loc[g] if g in hs.index else None)
    if isinstance(row, pd.DataFrame): row = row.iloc[0]
    nm = "" if row is None or pd.isna(row["name"]) else row["name"]
    ent = "" if row is None or pd.isna(row["entrez_id"]) else row["entrez_id"]
    loc = "" if row is None or pd.isna(row["location"]) else row["location"]
    lfc = float(de.lfc.get(g, np.nan)); padj = float(de.padj.get(g, np.nan))
    index.append([g, ens, nm, ent, loc, "AR" if g in AR else "NE" if g in NE else "",
                  round(float(pct_all[gi]) * 100, 1), round(float(mean_all[gi]), 3),
                  None if np.isnan(lfc) else round(lfc, 2), round(float(cor_ar[gi]), 3), round(float(cor_ne[gi]), 3)])
json.dump({"cols": ["symbol","ensembl","name","entrez","location","program","pct","mean","lfc_NE_vs_AR","r_AR","r_NE"], "rows": index, "shard": SHARD},
          open(OUT + "/genes.json", "w"), separators=(",", ":"))
# expression shards: uint8 dense, per-gene scale
for b in range(0, G, SHARD):
    blk = X[:, b:b+SHARD]; scl = np.maximum(blk.max(0), 1e-6) / 255
    q = np.round(blk / scl).astype(np.uint8).T.copy()
    open(f"{OUT}/expr/{b//SHARD:03d}.bin", "wb").write(q.tobytes())
    det = {}
    for j in range(blk.shape[1]):
        gi = b + j; g = genes[gi]
        det[g] = {"scale": float(scl[j]), "max": round(float(mx[gi]), 3),
                  "reg_mean": [round(float(v), 3) for v in rm[:, gi]], "reg_pct": [round(float(v)*100, 1) for v in rp[:, gi]],
                  "cl_mean": [round(float(v), 3) for v in cm[:, gi]], "cl_pct": [round(float(v)*100, 1) for v in cp[:, gi]],
                  "padj": None if pd.isna(de.padj.get(g, np.nan)) else float(f"{de.padj.get(g):.3g}"),
                  "cor": [[genes[k], round(float(v), 3)] for k, v in zip(topi[gi], topv[gi])]}
    json.dump(det, open(f"{OUT}/detail/{b//SHARD:03d}.json", "w"), separators=(",", ":"))
st = json.load(open("stats.json"))
sha = hashlib.sha256(open("analysis.py", "rb").read()).hexdigest()[:12]
try: commit = subprocess.check_output(["git", "-C", "portal", "rev-parse", "--short", "HEAD"]).decode().strip()
except Exception: commit = "n/a"
meta = {"n": int(ad.n_obs), "G": G, "img": {"w": img.width, "h": img.height}, "spot_px": sf["spot_diameter_fullres"] * s,
        "um_per_px": um_px / s, "x": np.round(xy[:, 0], 1).tolist(), "y": np.round(xy[:, 1], 1).tolist(),
        "barcode": ad.obs_names.tolist(), "region": [REG.index(v) for v in region], "regions": REG, "region_n": rn,
        "leiden": [CL.index(v) for v in cl], "clusters": CL, "cluster_n": cn,
        "AR_z": np.round(ad.obs.AR_z.values, 2).tolist(), "NE_z": np.round(ad.obs.NE_z.values, 2).tolist(),
        "provenance": {"geo": "GSE230282", "gsm": "GSM7211257 (EHU-W3)", "pmid": "37240308", "doi": "10.3390/ijms24108955",
                       "platform": "10x Visium CytAssist, FFPE, human prostate", "spots_in_tissue": st["spots_in_tissue"],
                       "spots_qc": st["spots_after_qc"], "genes": G, "median_umi": st["median_umi"], "median_genes": st["median_genes"],
                       "analysis_sha256": sha, "repo_commit": commit, "scanpy": sc.__version__, "z_threshold": st["z_threshold"],
                       "AR_genes": st["AR_genes_used"], "NE_genes": st["NE_genes_used"]}}
json.dump(meta, open(OUT + "/meta.json", "w"), separators=(",", ":"))
print("genes", G, "regions", dict(zip(REG, rn)))
