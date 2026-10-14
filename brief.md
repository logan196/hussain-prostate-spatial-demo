# Hussain call brief: AR → NE lineage in space (Wed Oct 14, 2026, 11:00 ET)

Live portal: see README for the URL · Script: `analysis.py` · Numbers: `stats.json` · Figures: `figures/`

## 1. Dataset (one, public, verified downloadable)
**GEO GSE230282 / GSM7211257 (EHU-W3).** 10x Visium CytAssist, FFPE prostate tissue from a 78-year-old man with **de novo metastatic NEPC co-existing with AR-positive adenocarcinoma (ARPC)** in the same section.
- Paper: Watanabe et al., *IJMS* 2023, [PMID 37240308](https://pubmed.ncbi.nlm.nih.gov/37240308/).
- I downloaded the matrix, positions, scalefactors and H&E directly from the GEO FTP (131 MB matrix).
- **Why this one:** it is the only public, directly downloadable human spatial transcriptomics section I found with large AR-high and NE-high tumor territories side by side. That lets us measure lineage segregation, overlap, and boundaries in one tissue.
- **Not a perfect fit:** this is de novo NEPC, not treatment-emergent, and there is no treatment history.
- **SHAPE cohort (Beltran lab):** verified as conference abstracts only (AACR 2026 PR027 and PCF 2025/2026). It covers 25 Visium samples across adeno, mixed CRPC/NEPC and NEPC, with treatment history. No public data release found. It is the dataset that would actually fit her question.

## 2. Her key claims relevant here
From her *Cancer Res Commun* 2026 paper (corresponding author); details in `hussain_profile.md`:
1. Prior ADT+ARPI (HR 3.04) or chemotherapy (HR 1.75) is linked to a **shorter time to t-SCNC**.
2. TP53 alterations are already present in 74% of baseline adenocarcinomas.
3. **RB1 loss is acquired** at transformation (10% → 76% in 21 paired biopsies).
4. Paired data show clonal continuity, i.e. the NE tumor descends from the adenocarcinoma.

## 3. What the spatial data shows
All numbers are computed by `analysis.py`.
- **QC:** 4,394 spots and 18,939 genes after QC. Median 23.3k UMI and 7.4k genes per spot.
- **Program scores:** AR program (11 genes: AR, KLK3, KLK2, TMPRSS2, NKX3-1, FKBP5, …) and NE program (14 genes: CHGA, CHGB, SYP, ASCL1, INSM1, NCAM1, ENO2, SOX2, …), scored with scanpy. "High" means z > 1.
- **Program calls:** 628 AR-high spots (14.3%), 890 NE-high (20.3%), and **22 dual AR+/NE+ (0.5% of all spots; 1.4% of program-high spots)**.
- **The two programs are mutually exclusive across the tissue.** AR and NE scores correlate at r = −0.43. Neighborhood enrichment (squidpy, 1,000 permutations) gives z = −30.9 for AR-high next to NE-high, versus +51.5 for AR next to AR and +79.6 for NE next to NE.
- **They occupy separate territories:**
  - Median distance from an NE-high spot to the nearest AR-high spot is **1.11 mm** (IQR 0.82–1.52 mm; max 2.44 mm).
  - Only 3.7% of NE spots are within 250 µm of an AR-high spot.
  - Only 4 of 628 pure AR-high spots directly touch an NE-high spot.
- **Gradient:** mean NE z rises steadily with distance from AR tissue: −0.55 at 0–100 µm, 0.27 at 0.5–1 mm, 1.19 at 1–2 mm. AR z falls in step.
- **Intermediate states:**
  - The 22 dual spots form one small cluster at the AR edge (median 79 µm from AR-high). Their neighborhood self-enrichment is z = 54.
  - They are **KLK3-high and CHGA-high but SYP-, ASCL1- and MKI67-low**. That looks like adenocarcinoma with neuroendocrine differentiation, not small-cell transition.
  - The NE-high region by contrast is proliferative (mean log MKI67 1.43 vs 0.10 in AR-high).
  - Treat these as a candidate, not a proven transition zone. Verify on H&E: the cluster sits near a dark pigment/edge artifact.
- **Unsupervised check:** Leiden clusters (12) recover the split. Clusters 0 and 4 are AR-dominant. Clusters 3 and 5 are NE-dominant (NE z 1.92 and 1.35).

**Bottom line for her hypothesis:**
- Spatially, NE and AR states in this tumor are distinct, sharply bounded territories, not a diffuse mixture. That fits a clonal "switch" model where a lineage-converted subclone expands, rather than widespread hybrid cells.
- It also means a single needle biopsy could easily sample only one state. That supports her case for paired and re-biopsy designs.
- These data **cannot** test treatment-sequence effects: one patient, de novo NEPC, no treatment annotation, no genomics.

## 4. Caveats to say up front
- n = 1 patient, one section, de novo NEPC (not t-SCNC), and no ADT/ARPI/chemo history.
- No DNA data here, so her RB1/TP53/PTEN claims are not testable. RB1 mRNA is not reduced in NE spots, consistent with the authors' report, but mRNA is not the same as genomic loss.
- 55-µm Visium spots contain several cells, so "dual" spots may be mixtures.
- Thresholds (z > 1) and gene sets are choices. The script exposes `--z` so they can be re-run live.

## 5. Questions to ask her
1. Her 21 paired biopsies: is any FFPE tissue left from baseline/t-SCNC pairs, especially mixed-histology cases? Spatial profiling (Visium HD/Xenium) of those pairs could map where RB1-loss subclones and NE programs arise, and whether ARPI exposure leaves a spatial "pre-NE" niche.
2. In her cohort, did mixed adeno/small-cell histology at transformation behave differently from pure t-SCNC? Would she want a spatial "AR–NE distance/mixing" metric as a candidate biomarker, e.g. for a re-biopsy trigger when PSA is discordant?
3. Would a prospective triplet-era cohort (ARASENS-like) with on-treatment biopsies be the right place to embed spatial profiling? Which endpoint matters most to her: time to t-SCNC, or response to platinum?
