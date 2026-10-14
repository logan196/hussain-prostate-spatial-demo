# Plan: Spatial map of treatment-driven lineage plasticity in prostate cancer (Hussain × Novaflow)
*Draft for the Oct 14, 2026 call. Prepared by Logan Beharry, founding engineer, Novaflow (YC S25).*

## 1. Question and hypothesis
**Question.** Does the order and intensity of prior therapy (ADT, then ARPI, then chemo/triplet) shape *where* and *how* adenocarcinoma turns into t-SCNC inside a tumor, and can that be read from tissue before the clinical switch?

**Starting point.** Dr. Hussain's six-center cohort of 140 cases [1] found:
- ADT+ARPI (HR 3.04) and chemotherapy (HR 1.75) shorten the time to t-SCNC.
- TP53 alterations are already present in 74% of baseline tumors.
- RB1 loss is acquired at transformation (10% → 76% in 21 paired biopsies).

Mechanistic work shows NE cells arise from KRT8+ luminal cells and then expand as clusters that depend on ASCL1 and require the in-vivo microenvironment [2]. FOXA2/NKX2-1 chromatin rewiring drives the switch [3]. CRPC also contains amphicrine (AR+/NE+) stem/progenitor states [4].

**Hypothesis.** ARPI exposure in a TP53-altered background selects spatially focal "pre-NE" niches at the AR boundary. In these niches, RB1-loss subclones pass through an amphicrine/ASCL1+ intermediate state and then expand as segregated NE territories. The size and distance of these niches measured at baseline predict time to t-SCNC.

**Pilot evidence (our demo, public GSE230282 [5]):**
- AR and NE programs sit in sharply segregated territories. Neighborhood enrichment for AR next to NE is z = −30.9, and the median NE-to-AR distance is 1.11 mm.
- Only 0.5% of spots carry both programs.
- The NE region is proliferative (MKI67).

**Caveat:** this is n = 1 de novo NEPC with no treatment history. It motivates the study but does not test the hypothesis.

## 2. What we build and what data we use
- **Cohort (core).** Her paired baseline/transformation FFPE biopsies, with treatment-sequence annotation [1]. Prioritize mixed adeno/SCNC cases.
  - Start with the 21 genomic pairs, then expand across her six centers.
  - Target 40–60 patients, roughly 100–150 sections.
- **Assays, per block:**
  - Visium HD (whole transcriptome, about 2 µm bins) for discovery and spatial CNV.
  - Xenium Prime 5K or CosMx 6K on serial sections for single-cell AR/NE/amphicrine calls.
  - CODEX/PhenoCycler protein panel for AR, PSA, INSM1, ASCL1, RB1, p53, Ki-67, and immune markers.
  - Existing clinical NGS for ground-truth TP53/RB1/PTEN.
- **Public anchors:**
  - GSE230282 [5].
  - Salachan et al. Visium CRPC/NEPC primaries [6].
  - SHAPE (Beltran lab: 25 Visium samples with treatment history; abstract only [7]). Approach the Beltran lab for access or collaboration.
- **Product.** A Novaflow "Lineage-Plasticity Atlas":
  - A static plus hosted portal, extending the current demo, with per-patient pre/post side-by-side, AR–NE axis overlays, clone maps, and treatment timelines.
  - Underneath it, an agentic pipeline that runs QC, annotation, CNV, and stats, versioned so every figure can be regenerated.

## 3. Core analyses
1. **Lineage scoring.** AR, NE, amphicrine, and double-negative scores, plus a continuous AR→NE pseudotime. Call transition states (KRT8+/ASCL1+, FOXA2/NKX2-1+).
2. **Niches.** Cellular neighborhoods (squidpy/CellCharter) and boundary-distance metrics. Stromal and immune composition at AR/NE interfaces.
3. **Clonality.** Spatial CNV and allele-specific CNV with CalicoST [8], validated against clinical NGS and following Erickson et al. [9]. Test whether NE territories descend from AR clones and where RB1/PTEN loss maps relative to them.
4. **Link to therapy.** Mixed models of niche metrics against prior regimen (ADT alone vs ADT+ARPI vs chemo vs triplet) and exposure duration, adjusted for site and center.
5. **Outcome.** Baseline spatial features against time to t-SCNC and OS (Cox models; internal cross-validation, with a held-out center for validation).

## 4. Publication plan
- **Target.** *Nature Medicine* (clinical cohort plus a predictive biomarker), with *Cancer Discovery* as the fallback. If the mechanism (clonal path, niches) dominates, aim for *Nature Cancer* or *Cancer Cell*. This matches where related work appears [2, 3].
- **Key figures:**
  - Fig 1: cohort, treatment timelines, and assay design.
  - Fig 2: paired pre/post spatial lineage maps.
  - Fig 3: transition and amphicrine niches at AR boundaries.
  - Fig 4: spatial clone phylogeography with RB1/TP53.
  - Fig 5: therapy sequence linked to niche metrics.
  - Fig 6: a baseline spatial score that predicts time to t-SCNC, plus the portal and a resource figure.
- **Timeline:**
  - Months 0–2: IRB/DUA, block pull, pilot of 6 pairs.
  - Months 2–6: profiling. The pipeline and portal run as data lands.
  - Months 6–9: analysis and validation.
  - Months 9–11: writing; preprint plus a public portal.
- **Roles:**
  - Hussain: clinical lead and senior/corresponding author (cohort, clinical questions, consortium).
  - Northwestern GU pathology: histology review and mixed-case calls.
  - Novaflow (Logan and team): computation, pipeline, portal. Co-first author for the computational lead.
  - Co-authors from the six centers.

## 5. What Novaflow adds
- **Speed.** Agentic QC, annotation, and analysis turn each batch around in days, not months.
- **Reproducibility.** Every number traces to versioned code (as in `analysis.py`/`stats.json`).
- **Portal.** Clinicians can explore cases with no code, and it becomes the paper's resource.
- **Scale.** Multi-platform harmonization (HD/Xenium/CosMx/CODEX) and reuse in future trial correlatives.

## 6. Asks for the call
1. Can we access FFPE from the 21 paired cases (and other mixed-histology cases), with treatment timelines and NGS? Who owns the blocks and data agreements across the six centers?
2. A pilot of 6 pairs (Visium HD plus Xenium on serial sections). Is there funding or a core to use (Lurie, PCF, or a DoD PCRP mechanism she knows well)? Should we co-write a grant?
3. Which endpoint matters most to her: time to t-SCNC, platinum response, or a re-biopsy trigger when PSA is discordant? And who is her pathology partner for histology calls?
4. Should we contact the Beltran lab (SHAPE) for access or collaboration?

## References (all verified)
1. Hussain M (corresp.) et al. Clinical and molecular profiles of treatment-emergent small cell neuroendocrine prostate cancers in the era of intensified systemic therapies. *Cancer Res Commun* 2026. doi:10.1158/2767-9764.CRC-26-0370
2. Romero R et al. The neuroendocrine transition in prostate cancer is dynamic and dependent on ASCL1. *Nat Cancer* 2024. doi:10.1038/s43018-024-00838-6
3. Lu X et al. NKX2-1 drives neuroendocrine transdifferentiation of prostate cancer via epigenetic and 3D chromatin remodeling. *Nat Genet* 2025. doi:10.1038/s41588-025-02265-4
4. Beshiri ML et al. Stem cell dynamics and cellular heterogeneity across lineage subtypes of castrate-resistant prostate cancer. *Stem Cells* 2024. doi:10.1093/stmcls/sxae025, PMID 38563224
5. Watanabe R et al. Spatial gene expression analysis reveals characteristic gene expression patterns of de novo NEPC coexisting with ARPC. *Int J Mol Sci* 2023. doi:10.3390/ijms24108955, PMID 37240308 (GEO GSE230282)
6. Salachan PV et al. Spatial whole transcriptome profiling of primary tumor from patients with metastatic prostate cancer. *Int J Cancer* 2023.
7. Bakht MK, …, Beltran H. Spatial heterogeneity atlas of prostate cancer evolution (SHAPE). *Cancer Res* 2026;86(2 Suppl):PR027 (abstract).
8. Ma C et al. Inferring allele-specific copy number aberrations and tumor phylogeography from spatially resolved transcriptomics (CalicoST). *Nat Methods* 2024;21:2239–2247. doi:10.1038/s41592-024-02438-9, PMID 39478176
9. Erickson A et al. Clonal phylogenies inferred from bulk, single cell, and spatial transcriptomic analysis of epithelial cancers. *PLoS One* 2025. doi:10.1371/journal.pone.0316475
