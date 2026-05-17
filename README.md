# VEN Fatigue Hypothesis: Social Cognition Paradox in Long-Duration Spaceflight

**Paper:** "The Social Cognition Paradox in Long-Duration Spaceflight: A VEN Fatigue Hypothesis for Duration-Dependent Emotion Recognition Decline"  
**Author:** Esila Keskin, University of the West of England, Bristol  
**Preprint:** _arXiv / bioRxiv link here once submitted_

---

## Overview

This repository contains all analysis code and data (where distributable) for the paper.
The study derives and tests the **VEN Fatigue Hypothesis**: that Von Economo neurons (VENs),
optimised for fast social decisions, are overactivated by gravity-disrupted social cues during
spaceflight, producing a duration-dependent threshold beyond which myelination compensation is
insufficient and Emotion Recognition Task (ERT) speed selectively declines.

Four independent datasets are tested:
- Cognitive performance: NASA Twins Study (340-day, N=1) and Dev et al. 2024 (6-month ISS, N=24)
- ISS rodent frontal cortex transcriptomics: GSE239336 (GeoMx DSP)
- Ground-based spaceflight analogue: OSD-202 (NASA OSDR)
- Human iPSC-derived cortical organoids on ISS: GSE259421

---

## Repository structure

```
ven-spaceflight-cognition/
  main.tex                       -- Paper LaTeX source
  references.bib                 -- Bibliography
  figures/                       -- Publication figures (PNG/PDF)
    fig1_ert_paradox.png
    fig2_domain_specificity.png
    fig3_molecular_dissociation.png
    fig4_organoid_permutation.png
  data/raw/
    twins_cognitive_heatmap.csv  -- Digitised from Garrett-Bakelman et al. 2019, Fig 10B
    dev2024_raw_scores.csv       -- From Dev et al. 2024, Table 2
    GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt  -- ISS rodent DE analysis (NCBI GEO)
    GLDS-202_...csv              -- OSD-202 (NASA OSDR) [gitignored, 134 MB]
  results/                       -- Analysis outputs (JSON, CSV)
    gse239336_ven_signature.json
    osd202_ven_signature.json
    VEN_panel_Combined.csv
    VEN_panel_NoMicroglia.csv
    VEN_panel_WithMicroglia.csv
    ...
  step1_create_cognitive_csvs.py -- Creates cognitive data CSVs
  step2_run_analysis.py          -- Runs ISS rodent and OSD-202 analyses
  step2_v2_run_analysis.py       -- Runs human organoid analysis (GSE259421)
  make_paper_figures.py          -- Generates all four publication figures
  ven_organoid_analysis.py       -- Supporting organoid analysis utilities
```

---

## Data availability

| Dataset | Source | Size | Status |
|---------|--------|------|--------|
| `twins_cognitive_heatmap.csv` | Digitised from Garrett-Bakelman et al. 2019, Science | 1 KB | Included |
| `dev2024_raw_scores.csv` | Dev et al. 2024, npj Microgravity, Table 2 | 1 KB | Included |
| `GSE239336_FCT_GCvsFLT-SAL_DEanalysis.txt` | NCBI GEO: GSE239336 | ~1 MB | Included |
| `GLDS-202_rna_seq_...csv` | NASA OSDR: OSD-202 | 134 MB | Gitignored -- download from [NASA OSDR](https://osdr.nasa.gov) |
| `GSE259421_all_counts.txt` | NCBI GEO: GSE259421 | 45 MB | Gitignored -- download from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259421) |

---

## Reproducing the analysis

### 1. Install dependencies

```bash
pip install numpy pandas scipy matplotlib
```

### 2. Download large data files (gitignored)

Download `GLDS-202_rna_seq_differential_expression_GLbulkRNAseq.csv` from NASA OSDR (OSD-202)
and `GSE259421_all_counts.txt` from NCBI GEO and place them at the paths above.

### 3. Run analysis pipeline

```bash
# Step 1: create cognitive performance CSVs from raw values
python step1_create_cognitive_csvs.py

# Step 2a: ISS rodent and OSD-202 analyses
python step2_run_analysis.py

# Step 2b: Human organoid analysis (GSE259421)
python step2_v2_run_analysis.py

# Step 3: Generate all four publication figures
python make_paper_figures.py
```

Figures are saved to `figures/` as both `.png` (150 dpi) and `.pdf` (vector).

---

## VEN gene panel

The 31-gene panel was defined in [Keskin 2026 (arXiv:2604.09229)](https://arxiv.org/abs/2604.09229)
prior to accessing any of the four datasets:

| Category | Genes |
|----------|-------|
| Myelination (n=7) | MBP, MOG, PLP1, MAG, CNP, MOBP, ERMN |
| Fast Signalling (n=7) | SCN1A, KCNQ2, ANK3, NEFH, NEFM, NEFL, SNCG |
| Social Circuit (n=6) | OXTR, AVPR1A, HTR2A, DRD1, CHRM1, GABRB2 |
| Layer V Projection (n=5) | FEZF2, BCL11B, TBR1, SATB2, CUX1 |
| Metabolic Support (n=6) | VDAC1, ATP2B2, SLC17A7, SNAP25, SYP, NRXN1 |

NOS1 (direct VEN biochemical marker) tracked individually outside the panel.

---

## Citation

If you use this code or data, please cite:

```
Keskin, E. (2026). The Social Cognition Paradox in Long-Duration Spaceflight:
A VEN Fatigue Hypothesis for Duration-Dependent Emotion Recognition Decline.
[Preprint]
```

---

## License

Code: MIT License. See LICENSE file.  
Data: subject to original dataset terms (NCBI GEO / NASA OSDR open data policies).
