# 🌠🔭 Assessing the Predictability of 𝛿 Scuti Variable Stars for Spacecraft Navigation

We assess the predictability of a large dataset of δ Scuti variable stars in detail. Few methods currently exist that can reliably determine the stability or predict the luminosity behavior of δ Scuti variable stars. To address this, we develop a computational framework to identify δ Scuti variable stars whose luminosity behavior can be accurately modeled as superpositions of sinusoidal functions. These predictive models provide a foundation for identifying δ Scuti variable stars suitable for use in practical applications such as autonomous deep-space navigation systems for spacecraft, which rely on such stars as navigational aid. This project provides the computational framework for **modeling the light curves of δ Scuti variable stars** outlined in the study. 

[![DOI](https://zenodo.org/badge/925039502.svg)](https://doi.org/10.5281/zenodo.16427747)

---

## 📂 Contents

- ``` bash
  final_code_base.py
  ```
  This is the **main script** that contains all the **methods** used to model, assess and analyze the light curves of **δ Scuti variable stars**. It includes:

  - Algorithms to model light curves as superpositions of sinusodial functions
  - R²_LSP and epsilon time error analysis
  - Detrending of TESS light curves and model accuracy assessment 
  - Plotting and interactive tools for light curves, target pixel files and periodograms


- ``` bash
  Master_Data_Sets_FULL
  ```
  This directory contains the **entire data set of analyzed δ Scuti variable stars**, including:

  - Cleaned light curve datasets
  - Star metadata (e.g., KIC IDs, TIC IDs, pulsation modes, magnitudes in Kepler photometric band, etc.)
  - Outputs (e.g., fitted parameters, NRMSE, epsilon, etc.)


- ``` bash
  Legacy_supporting_files
  ```
  Contains **legacy and non-essential scripts**, including:

  - Initial versions of the modeling code
  - Alternate modeling approaches
  - Additional plotting and visualization tools

  ⚠️ **Files within** `Legacy_supporting_files` **may not be stable**

---

## 📧 Points of Contact

- **Corresponding Author:** Ahmed Khan  ahmedk2@illinois.edu
- **Institution:** Department of Aerospace Engineering, University of Illinois at Urbana-Champaign
- **Collaborators:** Tiger Hou  linyhi2@illinois.edu, Siegfried Eggl  eggl@illinois.edu

---

## 💻 Usage

- ✅ **Actively Maintained** – The package is current and stable for general use.
- 🚫  **Not a Standalone Package** – To run the scripts, clone the repository and run them directly within Python. We recommend using an Anaconda environment. 
- ⚠️ **Requires Manual Data Handling** – Inputs will need to be directly inputed by the user.

---

## 🚀 Quickstart

### 🔧 Prerequisites

Make sure the following are installed:

```bash
pip install numpy pandas matplotlib scipy lightkurve scipy astropy unpopular scienceplots astroquery 
```

### ▶️ Running the Code

From the root directory:

**Run:**
```bash
python final_code_base.py
```

### 🌌 Sample Usage

Here is an example of using the framework within `final_code_base.py` to model a δ Scuti variable star:

```python
import re
import numpy as np
import matplotlib.pyplot as pt 
import lightkurve as lk 
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, lfilter
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from astroquery.mast import Tesscut
from astroquery.mast.utils import parse_input_location
import unpopular
import scienceplots

# Corresponds to KIC 3429637
name_star = '3429637'

# Generate the predictive model
predictive_model_array, light_curve, predictive_model_string = getCompositeSine2_deep(name_star) 

# Plot the composite function, residuals and light curve
plotsidebyside_deep(name_star)

# Extract epsilon values
epsilon_array, standard_deviation, slope_fit, P_MAX = get_epsilon_value(name_star, predictive_model_string)

# Produce the R^2_LSP and spectral residual value between the light curve and the predictive model

spectral_residuals, R2_LSP = SpectralResiduals(name_star, predictive_model_string)

```

---

## 📜 Background Documentation & Links

- [Gerald Handler; Delta Scuti Variables. AIP Conf. Proc. 16 September 2009; 1170 (1): 403–409.] (https://doi.org/10.1063/1.3246528)
- [Hou, L., Bansal, I., Davis, C., & Eggl, S. (2025). Position and Time Determination without Prior State Knowledge via Onboard Optical Observations of Delta Scuti Variable Stars. arXiv] (https://doi.org/10.48550/arXiv.2406.17609)

---

## 📁 Summary of Files: 

| Files                      | Purpose                        |
| -------------------------- | ------------------------------ |
| `final_code_base.py`       | Core Modeling Framework        |
| `Master_Data_Sets_FULL/`   | Complete Datasets              |
| `Legacy_supporting_files/` | Non-Essential, Legacy Scripts  |



