"""Independent recomputation of the v2 in-browser hypothesis tests.

Reads the same served files the browser reads (atlas/meta.json, vitessce/spots.json, vitessce/data.zarr/X)
and the result records exported from the browser (window.__results), then recomputes every number:
  - squidpy 1.8.2: spatial_neighbors(radius) + nhood_enrichment, co_occurrence
  - scipy: mannwhitneyu, cKDTree nearest distances
  - an exact replica of the browser's seeded permutation stream (mulberry32 + Fisher-Yates), so permutation
    p-values and null moments can be compared exactly, not only within Monte Carlo error.
Usage: python verify/verify_tests.py <repo_dir> <browser_results.json> <out_prefix>
"""
import json, sys, math
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

R, RES, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MT = json.load(open(f"{R}/atlas/meta.json")); SP = json.load(open(f"{R}/vitessce/spots.json"))
N = len(SP["barcode"]); x = np.array(MT["x"], float); y = np.array(MT["y"], float)
d1, _ = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2); pitch = float(np.median(d1[:, 1])); UM = 100.0 / pitch
X = x * UM; Y = y * UM; XY = np.c_[X, Y]

def mulberry32(seed):
    a = seed & 0xFFFFFFFF
    def nxt():
        nonlocal a
        a = (a + 0x6D2B79F5) & 0xFFFFFFFF; t = a
        t = ((t ^ (t >> 15)) * (t | 1)) & 0xFFFFFFFF
        t = (t ^ ((t + (((t ^ (t >> 7)) * (t | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0
    return nxt
def perms(labels, n, seed):
    rnd = mulberry32(seed); lab = labels.copy()
    for _ in range(n):
        for i in range(len(lab) - 1, 0, -1):
            j = int(rnd() * (i + 1)); lab[i], lab[j] = lab[j], lab[i]
        yield lab

def zarr_gene(sym):
    import zarr
    z = zarr.open(f"{R}/vitessce/data.zarr", mode="r"); var = [str(s) for s in z["var/_index"][:]]
    return np.asarray(z["X"][:, var.index(sym)], dtype=np.float64)

def group_mask(spec, other=None):
    k = spec["kind"]; op = spec.get("op", ">"); th = float(spec.get("threshold", 0) or 0)
    cmp = {">": np.greater, ">=": np.greater_equal, "<": np.less, "<=": np.less_equal}
    st = np.array(SP["state"]); reg = np.array(SP["region"]); lei = np.array([str(v) for v in MT["leiden"]])
    if k == "program":
        v = spec["value"].lower().replace(" ", "")
        m = reg == 0 if v.startswith("ar") else reg == 5 if v.startswith("ne") else reg == 2 if "dual" in v else ~np.isin(reg, [0, 5, 2])
    elif k == "state":
        v = spec["value"].lower()
        c = -1 if "abstain" in v else 2 if ("amphi" in v or "dual" in v) else 3 if "strom" in v else 1 if (v.startswith("ne") or "neuro" in v) else 0 if (v.startswith("ar") or "adeno" in v) else 4
        m = st == c
    elif k == "region": m = reg == [r.lower() for r in MT["regions"]].index(spec["value"].lower())
    elif k == "leiden": m = lei == "".join(ch for ch in str(spec["value"]) if ch.isdigit())
    elif k == "gene": m = cmp[op](zarr_gene(spec["gene"]), th)
    elif k == "score": m = cmp[op](np.array(MT["NE_z"] if str(spec.get("score", "")).upper().startswith("N") else MT["AR_z"]), th)
    elif k == "rest": m = ~other
    else: raise ValueError(k)
    for s in spec.get("and", []) or []: m = m & group_mask(s, other)
    return m

def labels_for(ma, mb):
    lab = np.full(N, -1, np.int32); lab[ma & ~mb] = 0; lab[mb & ~ma] = 1
    cats = ["A", "B"]; rest = (~ma) & (~mb); both = ma & mb
    if rest.any(): lab[rest] = len(cats); cats.append("rest")
    if both.any(): lab[both] = len(cats); cats.append("A and B")
    return lab, cats

def sq():
    sys.path.insert(0, "/tmp/sqt"); import sqstub  # noqa: light import shim, squidpy itself is unmodified
    import squidpy, anndata
    sys.meta_path[:] = [f for f in sys.meta_path if type(f).__module__ != "sqstub"]  # stop shimming once squidpy is loaded
    return squidpy, anndata

rows = []; recs = json.load(open(RES))
for rec in recs:
    t = rec["test"]; js = rec["result"]; P = rec["params"]
    ma = np.zeros(N, bool); ma[rec["idxA"]] = True; mb = np.zeros(N, bool); mb[rec["idxB"]] = True
    ga = group_mask(rec["spec"]["group_a"]); gb = group_mask(rec["spec"]["group_b"], ga)
    grp_ok = bool((ga == ma).all() and (gb == mb).all())
    lab, cats = labels_for(ma, mb); k = len(cats)
    row = {"id": rec["id"], "test": t, "A": rec["a"]["label"], "B": rec["b"]["label"], "nA_used": int((lab == 0).sum()), "nB_used": int((lab == 1).sum()), "nA_ui": rec["a"].get("n_used"), "nB_ui": rec["b"].get("n_used"), "groups_match": grp_ok}
    if t == "wilcoxon":
        v = zarr_gene(rec["gene"]); a = v[lab == 0]; b = v[lab == 1]
        r = mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
        row.update(stat="U", ui=js["U1"], py=float(r.statistic), ui_p=js["p"], py_p=float(r.pvalue), engine="scipy.stats.mannwhitneyu",
                   ok=abs(js["U1"] - r.statistic) < 1e-6 and (abs(js["p"] - r.pvalue) <= 1e-9 * max(r.pvalue, 1e-300) + 1e-300 or abs(math.log10(max(js["p"], 1e-320)) - math.log10(max(r.pvalue, 1e-320))) < 1e-6))
        row["auc_ui"] = js["auc"]; row["auc_py"] = float(r.statistic / (len(a) * len(b)))
    elif t == "nhood_enrichment":
        squidpy, anndata = sq()
        ad = anndata.AnnData(X=np.zeros((N, 1), np.float32)); ad.obsm["spatial"] = XY.copy()
        ad.obs["grp"] = pd.Categorical([cats[i] for i in lab], categories=cats)
        squidpy.gr.spatial_neighbors(ad, coord_type="generic", radius=float(P["radius_um"]))
        nE = int(ad.obsp["spatial_connectivities"].nnz)
        res = squidpy.gr.nhood_enrichment(ad, cluster_key="grp", n_perms=int(P["n_perms"]), seed=0, copy=True, show_progress_bar=False)
        z_sq = float(res.zscore[0, 1]); cnt_sq = float(res.counts[0, 1])
        # exact replica of the browser null
        A = ad.obsp["spatial_connectivities"].tocsr(); rr, cc = A.nonzero()
        def cnt(l): m = np.zeros((k, k)); np.add.at(m, (l[rr], l[cc]), 1); return m[0, 1]
        nul = np.array([cnt(p) for p in perms(lab, int(P["n_perms"]), 0)]); obs = cnt(lab)
        z_rep = (obs - nul.mean()) / nul.std(); pdep = (1 + (nul <= obs).sum()) / (len(nul) + 1); penr = (1 + (nul >= obs).sum()) / (len(nul) + 1)
        row.update(stat="A-B neighbour pairs / z", ui=f'{js["obs"]} / {js["zPair"]:.3f}', py=f"{cnt_sq:.0f} / {z_sq:.3f} (squidpy, own RNG)", replica=f"{obs:.0f} / {z_rep:.3f}",
                   ui_p=f'enr {js["pEnriched"]:.4g} dep {js["pDepleted"]:.4g}', py_p=f"enr {penr:.4g} dep {pdep:.4g} (replica)", edges_ui=js["edges"], edges_py=nE, engine="squidpy.gr.nhood_enrichment 1.8.2",
                   ok=(js["obs"] == cnt_sq) and js["edges"] == nE and abs(js["zPair"] - z_rep) < 1e-6 and abs(js["zPair"] - z_sq) < max(0.1 * abs(z_sq), 0.5) and abs(js["pDepleted"] - pdep) < 1e-12 and abs(js["pEnriched"] - penr) < 1e-12,
                   z_squidpy=z_sq, z_ui=js["zPair"], z_diff_vs_squidpy=js["zPair"] - z_sq)
    elif t == "proximity":
        scope = P.get("null_scope", "all")
        def stat(l):
            A = np.where(l == 0)[0]; B = np.where(l == 1)[0]
            return float(np.mean(cKDTree(XY[B]).query(XY[A])[0]))
        def nulls(n):
            if scope != "ab":
                for p in perms(lab, n, 0): yield p
            else:
                pos = np.where((lab == 0) | (lab == 1))[0]; l2 = lab.copy()
                for sub in perms(lab[pos].copy(), n, 0): l2[pos] = sub; yield l2
        obs = stat(lab); nul = np.array([stat(p) for p in nulls(int(P["n_perms"]))])
        pc = (1 + (nul <= obs).sum()) / (len(nul) + 1); pf = (1 + (nul >= obs).sum()) / (len(nul) + 1); z = (obs - nul.mean()) / nul.std()
        row.update(stat=f"mean nearest A-B (um) / z [null: {scope}]", ui=f'{js["obs"]:.3f} / {js["z"]:.3f}', py=f"{obs:.3f} / {z:.3f}", ui_p=f'closer {js["pCloser"]:.4g} farther {js["pFarther"]:.4g}', py_p=f"closer {pc:.4g} farther {pf:.4g}",
                   engine="scipy cKDTree + seeded replica", ok=abs(js["obs"] - obs) < 1e-6 and abs(js["z"] - z) < 1e-6 * max(1, abs(z)) and abs(js["pCloser"] - pc) < 1e-12 and abs(js["pFarther"] - pf) < 1e-12)
    elif t == "co_occurrence":
        squidpy, anndata = sq()
        ad = anndata.AnnData(X=np.zeros((N, 1), np.float32)); ad.obsm["spatial"] = XY.copy()
        ad.obs["grp"] = pd.Categorical([cats[i] for i in lab], categories=cats)
        occ, interval = squidpy.gr.co_occurrence(ad, cluster_key="grp", interval=np.array(P["interval_um"], float), copy=True)
        curve_sq = occ[0, 1, :]; dmax = float(np.max(np.abs(curve_sq - np.array(js["curve"]))))
        d = float(P["distance_um"]); tree = cKDTree(XY); pairs = tree.query_pairs(d, output_type="ndarray")
        I = np.r_[pairs[:, 0], pairs[:, 1]]; J = np.r_[pairs[:, 1], pairs[:, 0]]
        def st(l):
            c = np.zeros((k, k)); np.add.at(c, (l[I], l[J]), 1); row_ = c.sum(0); tot = row_.sum()
            return (c[1, 0] / row_[1]) / (row_[0] / tot)
        obs = st(lab); nul = np.array([st(p) for p in perms(lab, int(P["n_perms"]), 0)])
        ph = (1 + (nul >= obs).sum()) / (len(nul) + 1); pl = (1 + (nul <= obs).sum()) / (len(nul) + 1); z = (obs - nul.mean()) / nul.std()
        idx = list(P["interval_um"][1:]).index(d) if d in P["interval_um"][1:] else None
        row.update(stat=f"co-occurrence ratio @ {d:.0f} um / z", ui=f'{js["obs"]:.5f} / {js["z"]:.3f}', py=f"{obs:.5f} / {z:.3f}" + (f" (squidpy {curve_sq[idx]:.5f})" if idx is not None else ""),
                   ui_p=f'higher {js["pHigher"]:.4g} lower {js["pLower"]:.4g}', py_p=f"higher {ph:.4g} lower {pl:.4g}", curve_max_abs_diff_vs_squidpy=dmax, pairs_ui=js["pairs"], pairs_py=int(len(I)),
                   engine="squidpy.gr.co_occurrence 1.8.2 + seeded replica", ok=dmax < 1e-4 and abs(js["obs"] - obs) < 1e-9 and abs(js["pHigher"] - ph) < 1e-12 and js["pairs"] == len(I))
    rows.append(row); print(json.dumps(row, default=str))
df = pd.DataFrame(rows); df.to_csv(OUT + ".csv", index=False)
json.dump({"pitch_px": pitch, "um_per_px": UM, "v1_um_per_px": MT["um_per_px"], "rows": rows}, open(OUT + ".json", "w"), indent=1, default=str)
print("ALL OK" if all(r.get("ok") and r["groups_match"] and r["nA_ui"] in (None, r["nA_used"]) and r["nB_ui"] in (None, r["nB_used"]) for r in rows) else "MISMATCH")
