This repository accompanies the paper *"Drift Modeling for NDIR Gas Sensors: Separating Environmental Effects and Instrumental Drift"*. It contains code written in both R and Python to generate figures and analyze synthetic and real-world sensor data.

## Contents

- **Conceptual Illustration (Fig. 2)**  
  Jupyter notebook: `concept of drift/stochastic process.ipynb`

- **Synthetic Example (removed in revision)**  
  R script: `PoF experiments/synthetic example/simulation strain gauge.r`

- **In-field Example**  
  - The data used in this study are not publicly available because they are not owned by us and we do not have the rights to share them.
  - Data pre-processing: `in-field experiment/.with/outlier_remove_all.R`
  - GAM fit for each individual sensor and corresponding plots: `in-field experiment/GAM_for_all_new.R`
  - Improved baseline calibration: `in-field experiment/baseline T humidity.R`
  - Model comparison: `in-field experiment/.with/without interaction bam.R`
  - Drawing plots that combine the effects from all sensors: `in-field experiment/plot for all.r`

## Requirements

### Python (for the Jupyter notebook)
- Python 3
- Jupyter
- Required libraries:
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `scienceplots`

### R
- `mgcv`
- `ggplot2`, `dplyr`, `tidyr`

---
For questions or feedback, please contact: [yangch@kth.se]