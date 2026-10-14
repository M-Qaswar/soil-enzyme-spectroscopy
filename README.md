# Soil Enzyme Spectroscopy — Reproducible Analysis Package

Analysis code accompanying:

**Fusion Strategies for VIS–NIR and MIR Spectroscopy for Predicting Soil Microbial Enzyme Activities**
Qaswar, M., Kuramae, E.E., Mouazen, A.M. *Computers and Electronics in Agriculture* (under review).

Archived release: https://doi.org/10.5281/zenodo.23246399

## Files

- `RRSEL_pipeline.R` — PLSR preprocessing and model development, Kennard–Stone partitioning, Fusion-PCA, RRSEL stacking, VIP and validation bootstrap, performance metrics, paired statistical comparisons, Bland–Altman calculations.
- `RF_benchmark.py` — random-forest benchmark with calibration-only hyperparameter optimisation.

Both scripts are annotated. The methodological rationale, preprocessing choices, partitioning, hyperparameter grid and leakage safeguards are described in the manuscript (Sections 2.3–2.9 and Supplementary Table S1); this file documents only what is needed to run the code.

## Input data

CSV files must contain:

1. `sampleID`
2. `BG`, `ACP`, `ARS` — enzyme activities in µg PNP g-1 dry soil h-1
3. `ALP` — retained for compatibility with the supplied data format; not used
4. VIS–NIR reflectance variables from 305–1700 nm at 1-nm spacing
5. MIR absorbance variables from 4000–650 cm-1

Spectral columns are identified from their numeric wavelength/wavenumber names. If the original CSV uses another naming convention, rename those columns to numeric values before running the scripts.

Sample numbers: BG is modelled on 111 samples (78 calibration / 33 validation); ACP on 92 (67 / 25); ARS on 92 (69 / 23). The random forest benchmark uses the 92 samples common to all three enzymes, with its own partition of 69 / 23. The code does not create or duplicate observations.

## Software

R 4.5.1 with `prospectr`, `caret`, `pls`, `glmnet`, `zoo`, `ggplot2`, `patchwork`.
Python 3.11 with `numpy`, `pandas`, `scipy`, `scikit-learn`.

Random seeds are fixed in both scripts (principal seed 125).

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

## Data availability

Raw soil spectra and enzyme-activity data are not distributed with this package. The code permits methodological verification when the corresponding data are supplied; numerical identity with every published table value cannot be guaranteed from the code alone.
