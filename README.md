# Soil Enzyme Spectroscopy — Reproducible Analysis Package

This repository contains the two annotated source files implementing the analytical workflow described in the manuscript:

**Fusion Strategies for VIS–NIR and MIR Spectroscopy for Predicting Soil Microbial Enzyme Activities**

## Files

- `RRSEL_pipeline.R` — PLSR preprocessing/model-development pipeline, Kennard–Stone partitioning, Fusion-PCA, RRSEL stacking, VIP bootstrap, validation bootstrap, performance metrics, paired statistical comparisons, and Bland–Altman calculations.
- `RF_benchmark.py` — independently implemented random-forest benchmark with calibration-only five-fold hyperparameter optimisation.
- `README.md` — documentation and reproducibility notes.

## Data-design requirement

The revised manuscript reports 111 samples for BG and 92 samples for ACP and ARS, while Section 2.3 currently states that spectra were acquired for all 92 samples. These statements should be reconciled in the final manuscript.

The code does not create, duplicate, or otherwise invent observations to reconcile differences in sample counts. If 111 BG observations are used in the final analysis, the reported BG partition of 78 calibration and 33 validation observations can be reproduced. ACP and ARS use 92 observations. The RF benchmark uses the 92 observations common to BG, ACP, and ARS and therefore uses a 69/23 calibration–validation partition.

The origin of the additional 19 BG observations, if applicable to the final study design, should be explicitly described in the manuscript.

## Input data

CSV files must contain:

1. `sampleID`
2. `BG`
3. `ACP`
4. `ALP` (retained for compatibility with the supplied data format; not used)
5. `ARS`
6. VIS–NIR reflectance variables from 305–1700 nm at 1-nm spacing
7. MIR absorbance variables from 4000–650 cm-1

The code identifies spectral columns from their numeric wavelength/wavenumber names where possible. If the original CSV uses another naming convention, rename the spectral columns to numeric values before running the scripts.

Enzyme activities must be in:

`µg PNP g-1 dry soil h-1`

## Method implemented

### PLSR/RRSEL workflow

For each enzyme:

- VIS–NIR is converted from reflectance to apparent absorbance.
- Detector-junction correction is applied at 1024–1044 nm.
- The modelling ranges are VIS–NIR 394–1674 nm and enzyme-specific MIR ranges:
  - BG: 3954.7–853.6 cm-1
  - ACP: 3928.6–853.6 cm-1
  - ARS: 3954.7–853.6 cm-1
- The Table 2 preprocessing chain is selected using calibration-only LOOCV.
- Kennard–Stone is applied to VIS–NIR using the enzyme-specific PC counts and calibration fractions:
  - BG: 78/111 calibration, 33/111 validation; fraction 0.70; 11 PCs
  - ACP: 67/92 calibration, 25/92 validation; fraction 0.73; 10 PCs
  - ARS: 69/92 calibration, 23/92 validation; fraction 0.75; 11 PCs
- The resulting sample indices are applied unchanged to MIR.
- PLSR latent variables are searched from 1–20 by LOOCV within calibration only.
- M1 = VIS–NIR PLSR.
- M2 = MIR PLSR.
- M3 = direct VIS–NIR + MIR concatenation PLSR.
- M4 = independent calibration-only PCA of each spectral block, retaining the minimum number of PCs explaining at least 95% of calibration variance, followed by concatenation and PLSR.
- M5 = RRSEL using out-of-fold calibration predictions from M1–M3 and a ridge meta-learner. The meta-learner is trained only on calibration data.
- Validation samples are untouched until final prediction.
- VIP bootstrap uses 500 calibration resamples.
- Validation confidence intervals use 2,000 paired bootstrap resamples.
- Friedman and paired Wilcoxon tests use the same external-validation samples for M1–M5, with Holm correction.

### Random forest benchmark

The RF benchmark is separate from the PLSR/RRSEL pipeline.

- Only the 92 samples common to BG, ACP and ARS are used.
- A separate 75% Kennard–Stone partition is generated from VIS–NIR using the minimum number of PCs explaining at least 95% of VIS–NIR variance.
- This gives 69 calibration and 23 validation samples.
- VIS–NIR and MIR use the same RF sample indices.
- Preprocessing: VIS–NIR absorbance, Savitzky–Golay smoothing (window 5, polynomial 2), then first derivative (window 11, polynomial 2). MIR uses the corresponding common RF preprocessing specified in Section 2.4.
- Hyperparameters are tuned **only within the calibration set** using five-fold shuffled K-fold CV:
  - 500 trees
  - `max_features`: `sqrt(p)`, `0.1p`, `0.33p`
  - `max_depth`: unlimited, 10
  - `min_samples_leaf`: 1, 3, 5
- The selected model is refitted on all calibration samples and evaluated once on the 23 external validation samples.
- RF is not included in the formal five-model Friedman/Wilcoxon comparison.

## Leakage prevention

Operations that act independently on a single spectrum may be performed before partitioning because they estimate no parameters from other samples. Operations that estimate parameters across samples are calibration-only.

The implementation therefore keeps:

- centering/scaling calibration-only;
- PCA calibration-only;
- PLSR model selection calibration-only;
- RF hyperparameter optimisation calibration-only;
- ridge meta-model fitting calibration-only;
- external validation completely unavailable to model-development decisions.

The RRSEL meta-learner is fitted to leave-one-out out-of-fold predictions rather than in-sample base-model predictions.

## Reproducibility

The principal random seed is 125. The RF cross-validation uses a fixed seed of 125. The R workflow uses a separate resampling seed of 126 where applicable.

The code implements the workflow and permits methodological verification when the corresponding data are supplied. Because the original spectral and enzyme-activity data are not distributed with this package, numerical identity with every published table value cannot be guaranteed from the code alone.

## Software

R 4.5.1 with `prospectr`, `caret`, `pls`, `glmnet`, `zoo`, `ggplot2`, and `patchwork`.

Python 3.11 with `numpy`, `pandas`, `scipy`, and `scikit-learn`.

## Running

### R

```r
source("RRSEL_pipeline.R")
vis <- read.csv("visnir_spectra.csv", check.names = FALSE)
mir <- read.csv("mir_spectra.csv", check.names = FALSE)

res <- run_enzyme(vis, mir, "BG")
res$performance
res$vip
res$comparison
```

### Python

```python
import pandas as pd
from RF_benchmark import run_enzyme

vis = pd.read_csv("visnir_spectra.csv")
mir = pd.read_csv("mir_spectra.csv")

res = run_enzyme(vis, mir, "BG")
print(res["performance"])
print(res["best_params"])
```

The package intentionally does not contain raw soil spectra or enzyme-activity data.
