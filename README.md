# LBM-R
## Repository Overview

This repository contains all the code, data, and supporting materials related to the analysis published in the British Journal of Radiology (BJR) article titled: “Puri T, Blake GM.  Comparison of ten predictive equations for estimating lean body mass with dual-energy X-ray absorptiometry in older patients. British Journal of Radiology. 2022;95(1133):20210378"

### Contents

This repository includes 8 files and 1 folder, each serving a specific purpose:

1. `Article_BJR_20210378.pdf` : This is the full peer-reviewed scientific article published in BJR. It presents the results of analyses performed using the original DXA dataset (not included here) and the associated R scripts provided here.

2. `Article_BJR_20210378-Supplement.pdf` : This file contains the supplementary material that was originally published alongside the main article. As of now, this supplement is no longer accessible from the BJR website (Error 404), so a local copy is included for reference.

3. `LICENSE` : Specifies the licensing terms under GPL-3.0. Users are free to use, modify, and distribute this code under the terms of the GNU General Public License v3.

4. `README.md` : This file (you're reading it) provides a structured explanation of the repository contents and how to use them.

5. `simulate_synthetic_dxa_dataset.R` : R script used to generate synthetic DXA data for testing and validation. While the simulated data does not precisely replicate real clinical distributions, it is useful for verifying the logic and reproducibility of the analysis code.

6. `simulated_dxa_dataset.csv` : A mock DXA dataset with realistic-looking numerical patterns. This dataset does not represent real patient data but was randomly generated to mimic the structure of the original dataset (which was obtained from Erasmus Professor Glen Blake at King's College London, UK). The data in this file was used as input to `dxa_data_analysis.R` to generate the output figures stored in the  `results/` (folder).

7. `simulated_dxa_dataset_from_code.csv` : Synthetic dataset generated using the R script in `simulate_synthetic_dxa_dataset.R`. It serves as an example output of the simulation process.

8. `dxa_data_analysis.R` : This directory contains the core R script used for processing and modeling both real and synthetic DXA data, generating results comparable to those presented in the BJR article (`Article_BJR_20210378.pdf`) and its supplementary material (`Article_BJR_20210378-Supplement.pdf`).

9. `results/` (folder)
   A folder containing example result figures (in `.jpg` format) produced by the analysis code in `dxa_data_analysis.R` file. These images demonstrate typical outputs of the modeling and Bland-Altman comparisons implemented in the scripts. Correlation plots are printed to the console, while Bland–Altman plots are saved as image files.

### Notes

* The original DXA dataset is not publicly available due to patient confidentiality and data-sharing agreements. However, the structure and logic of the analysis can be replicated using the simulated data provided.
* This repository is provided in the spirit of transparency, reproducibility, and reuse, in accordance with the terms of the GNU General Public License v3 (GPL-3.0).






