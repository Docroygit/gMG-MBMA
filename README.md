# Regimen-stratified Bayesian model-based meta-analysis of targeted biologics in generalised myasthenia gravis

Code and extracted data for a Bayesian model-based meta-analysis (MBMA) of targeted biologics in generalised myasthenia gravis (gMG), published in *Naunyn-Schmiedeberg's Archives of Pharmacology*.

Nine placebo-controlled trials are analysed on two efficacy scales (MG-ADL and QMG). Dose-response is described by hierarchical nonlinear Emax models fitted with brms/Stan, with separate models for continuous and cyclical/weekly regimens and a placebo model per scale (six models in total).

| Trial | Drug |
|---|---|
| REGAIN | eculizumab |
| CHAMPION-MG | ravulizumab |
| RAISE | zilucoplan |
| MINT | inebilizumab |
| ADAPT | efgartigimod |
| VIVACITY-MG3 | nipocalimab |
| LUMINESCENCE | satralizumab |
| batoclimab RCT | batoclimab |
| MYCARING | rozanolixizumab |

## Repository layout

```
R/gMG_MBMA_analysis.R              full analysis script
data/gMG_MBMA_extraction.xlsx      extracted trial data
outputs/                           created on first run
```

## Data

`gMG_MBMA_extraction.xlsx` has three sheets:

- `MGADL_change` – MG-ADL change from baseline by trial, arm and time point
- `QMG_change` – QMG change from baseline
- `Tolerability_safety_profiles` – route, dosing schedule and safety/tolerability figures per drug, used in the benefit-risk module

## Requirements

R (>= 4.2) with a working Stan toolchain (Rtools on Windows), and the packages
`readxl`, `dplyr`, `brms`, `tidybayes`, `ggplot2`, `tidyr`, `stringr`, `patchwork`, `gridExtra`, `pracma`, `forcats` (`grid` ships with R).

```r
install.packages(c("readxl", "dplyr", "brms", "tidybayes", "ggplot2", "tidyr",
                   "stringr", "patchwork", "gridExtra", "pracma", "forcats"))
```

## Running

From the repository root:

```r
source("R/gMG_MBMA_analysis.R")
```

The script changes the working directory to `outputs/` and writes all figures, tables and model objects there. MCMC fitting is slow; expect several hours for the full run. Fitted models are saved after the QMG section (`all_fitted_models.RData`) and the whole workspace at the end (`gMG_MBMA_complete_workspace.RData`).

## Script structure

Sections are marked `# SECTION N:` in the script.

| Sections | Content |
|---|---|
| 1-13 | MG-ADL: data preparation, model fitting, diagnostics, VPCs, Emax, week-26 net benefit, time to MCID, AUEC |
| 14-21 | QMG: same pipeline |
| 22 | Save fitted models |
| 23-25 | Clinical trial simulation, treatment ranking (SUCRA), benefit-risk-convenience |
| 26-31 | Model comparison (LOO), prior sensitivity, leave-one-study-out, cross-scale concordance, mechanism contrast, heterogeneity and ED50 identifiability |
| 32-33 | Sensitivity analyses: MCID threshold, MINT placebo exclusion |
| 34-35 | Model summaries and final save |

## Model notes

- Every brms formula uses `se(se, sigma = TRUE)` so that reported standard errors enter the likelihood while residual variance is still estimated.
- Continuous regimens use cumulative dose as exposure; cyclical/weekly regimens use the on-period average dose.
- Random seeds are set in the script (`set.seed(2024)` and `seed = 123` in `brm()`).

## Citation

Please cite the accompanying article.
