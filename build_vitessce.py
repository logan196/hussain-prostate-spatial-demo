#!/usr/bin/env python3
"""Build Vitessce core data: full-gene AnnData zarr (v2), ensemble cell-state calls with abstention, hypothesis stats.
Run after analysis.py and build_atlas.py."""
import json, shutil, numpy as np, pandas as pd, scanpy as sc, anndata, scipy.sparse as sp, tifffile
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture
from PIL import Image
anndata.settings.zarr_write_format = 2; anndata.settings.allow_write_nullable_strings = False
OUT = "portal/vitessce"; ABSTAIN = 0.75
a = sc.read_h5ad("data/processed.h5ad"); m = json.load(open("portal/atlas/meta.json"))
X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X); X = X.astype(np.float32)
STR = [g for g in ["COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","MYH11","VIM","TAGLN","SPARC"] if g in a.var_names]
sc.tl.score_genes(a, STR, score_name="STR_score", random_state=0)
z = lambda v: (v - v.mean()) / v.std()
S = np.c_[z(a.obs.AR_score.values), z(a.obs.NE_score.values), z(a.obs.STR_score.values)]
LAB = ["AR adenocarcinoma", "NE / t-SCNC-like", "Amphicrine (AR+NE)", "Stroma", "Other"]
def rule(s, t=1.0, ta=0.5):
    ar, ne, st = s[:, 0], s[:, 1], s[:, 2]
    out = np.full(len(s), 4)
    out[(st > t) & (ar < t) & (ne < t)] = 3
    out[(ar > t) & (ne <= ta)] = 0; out[(ne > t) & (ar <= ta)] = 1
    out[(ar > ta) & (ne > ta) & ((ar > t) | (ne > t))] = 2
    return out
v1 = rule(S)                                              # voter 1: z-score thresholds
v2 = rule(S, t=0.75, ta=0.4)                              # voter 2: relaxed thresholds
A = a.obsp["spatial_connectivities"].astype(float); deg = np.asarray(A.sum(1)).ravel()
Sn = (A @ S + S) / (deg[:, None] + 1); v3 = rule(Sn)       # voter 3: spatially smoothed scores
gm = GaussianMixture(8, random_state=0, n_init=3).fit(S)  # voter 4: GMM components labelled by mean profile
v4 = rule(gm.means_, t=0.6, ta=0.3)[gm.predict(S)]
cl = a.obs.leiden.astype(str).values                       # voter 5: Leiden cluster mean profile
cm = {c: S[cl == c].mean(0) for c in np.unique(cl)}; v5 = np.array([rule(cm[c][None], t=0.6, ta=0.3)[0] for c in cl])
V = np.c_[v1, v2, v3, v4, v5]
cnt = np.stack([(V == k).sum(1) for k in range(5)], 1); maj = cnt.argmax(1); agree = cnt.max(1) / V.shape[1]
state = np.where(agree >= ABSTAIN, maj, -1)
names = np.array([LAB[k] if k >= 0 else "Abstained (low confidence)" for k in state])
a.obs["Cell state"] = pd.Categorical(names); a.obs["Program region"] = pd.Categorical([m["regions"][i] for i in m["region"]])
a.obs["Leiden"] = pd.Categorical(cl)
# hypothesis stats
d = a.obs.dist_to_AR_um.values; nonAR = m["region"] != np.array(0) if False else np.array(m["region"]) != 0
rho, p = spearmanr(d[nonAR], a.obs.NE_z.values[nonAR])
ne_lab = names == LAB[1]
hyp = {"spearman_dist_vs_NEz_nonAR": [round(float(rho), 3), float(f"{p:.2g}")], "n_nonAR": int(nonAR.sum()),
       "NE_state_dist_median_um": round(float(np.median(d[ne_lab])), 0),
       "NE_state_within_500um_of_AR": int((ne_lab & (d <= 500)).sum()), "NE_state_total": int(ne_lab.sum()),
       "amphicrine_n": int((names == LAB[2]).sum()), "amphicrine_dist_median_um": round(float(np.median(d[names == LAB[2]])), 0) if (names == LAB[2]).any() else None}
ens = {"labels": LAB, "voters": ["z>1 thresholds", "relaxed z>0.75", "spatially smoothed (6-NN)", "GMM (8 comp.)", "Leiden-cluster profile"],
       "abstain_below_agreement": ABSTAIN, "counts": {k: int(v) for k, v in pd.Series(names).value_counts().items()},
       "abstention_rate_pct": round(100 * float((state < 0).mean()), 2), "mean_agreement": round(float(agree.mean()), 3),
       "unanimous_pct": round(100 * float((agree == 1).mean()), 2), "stroma_genes": STR}
# matrix: all genes + derived features (prefixed)
extra = {"◆ state confidence": agree, "◆ AR score (z)": S[:, 0], "◆ NE score (z)": S[:, 1], "◆ stroma score (z)": S[:, 2], "◆ distance to AR (mm)": d / 1000}
Xf = np.c_[X, np.c_[[extra[k] for k in extra]].T].astype(np.float32); Xf = np.round(Xf, 2)
var = list(a.var_names) + list(extra)
b = anndata.AnnData(Xf, obs=pd.DataFrame(index=a.obs_names.astype(object)), var=pd.DataFrame(index=pd.Index(var, dtype=object)))
SHORT = {"AR adenocarcinoma": "AR adeno", "NE / t-SCNC-like": "NE / t-SCNC", "Amphicrine (AR+NE)": "Amphicrine", "Abstained (low confidence)": "Abstained"}
for c in ["Cell state", "Program region", "Leiden"]: b.obs[c] = pd.Categorical(a.obs[c].astype(str).map(SHORT).fillna(a.obs[c].astype(str)) if c == "Cell state" else a.obs[c].astype(str))
b.obsm["spatial"] = np.c_[m["x"], m["y"]].astype(np.float32); b.obsm["X_umap"] = a.obsm["X_umap"].astype(np.float32)
shutil.rmtree(OUT + "/data.zarr", ignore_errors=True); b.write_zarr(OUT + "/data.zarr", chunks=(b.n_obs, 40))
# H&E OME-TIFF (interleaved RGB)
im = np.array(Image.open("portal/atlas/tissue_hires.jpg").convert("RGB"))
# pyramidal, tiled, deflate-compressed: the viewer fetches only the tiles on screen
from PIL import Image as _I
cyx = im.transpose(2, 0, 1); levels = [cyx]
while max(levels[-1].shape[1:]) > 512:
    h, w = levels[-1].shape[1:]; levels.append(np.stack([np.array(_I.fromarray(c).resize((w // 2, h // 2), _I.LANCZOS)) for c in levels[-1]]))
with tifffile.TiffWriter(OUT + "/he.ome.tif", ome=True, bigtiff=False) as tw:
    opts = dict(photometric="minisblack", tile=(256, 256), compression="jpeg", compressionargs={"level": 75})
    tw.write(levels[0], subifds=len(levels) - 1, metadata={"axes": "CYX", "Channel": {"Name": ["R", "G", "B"]}}, **opts)
    for lv in levels[1:]: tw.write(lv, subfiletype=1, **opts)
# client-side per-spot data for live queries
json.dump({"barcode": list(a.obs_names), "state": [int(s) for s in state], "conf": np.round(agree, 2).tolist(), "dist": np.round(d).astype(int).tolist(),
           "region": m["region"], "regions": m["regions"], "labels": LAB, "ensemble": ens, "hyp": hyp},
          open(OUT + "/spots.json", "w"), separators=(",", ":"))
print(json.dumps(ens, indent=1)); print(hyp)
# re-write X with stronger compression (blosc zstd + bitshuffle), column chunks of 40 genes
import zarr, numcodecs, subprocess
g = zarr.open_group(OUT + "/data.zarr", mode="r+", zarr_format=2); attrs = dict(g["X"].attrs); del g["X"]
arr = g.create_array("X", shape=Xf.shape, chunks=(Xf.shape[0], 40), dtype="float32", fill_value=0.0,
                     compressor=numcodecs.Blosc(cname="zstd", clevel=9, shuffle=numcodecs.Blosc.BITSHUFFLE))
arr[:] = Xf; arr.attrs.update(attrs)
zarr.consolidate_metadata(OUT + "/data.zarr", zarr_format=2)
print(subprocess.check_output(["du", "-sh", OUT + "/data.zarr"]).decode())
