# 🌠🔭 Fourier-Based Predictive Modeling of δ Scuti Variable Stars Framework

We assess the predictability of a large dataset of δ scuti variable stars in detail. Few methods currently exist that can reliably determine the stability or predict the luminosity behavior of δ scuti variable stars. To address this, we develop a computational framework to identify δ scuti variable stars whose luminosity behavior can be accurately modeled as superpositions of sinusoidal functions. These predictive models provide a foundation for identifying δ scuti variable stars suitable for use in practical applications such as autonomous deep-space navigation systems for spacecraft, which rely on such stars as navigational aid. This project provides the computational framework for **modeling the light curves of delta scuti variable stars** outlined in the study. 


---

## 📚 Package Contents

- ``` bash
  final_code_base.py
  ```
  This is the **main script** that contains all the **methods** used to model, assess and analyze the light curves of **delta scuti variable stars**. It includes:

  - Fourier decomposition modeling
  - R_fft² and residual analysis
  - Preprocessing and smoothing routines
  - Plotting utilities for light curves

- ``` bash
  Legacy_supporting_files
  ```
  This directory contains the **entire set of analyzed observational data**, including:

  - Cleaned light curve datasets
  - Star metadata (KIC IDs, TIC IDs, pulsation modes, magnitude in  etc.)
  - Sample outputs (e.g., fitted parameters and residuals)

- ``` bash
  Master_Data_Sets_FULL
  ```
  Contains **archived and non-essential scripts**, including:

  - Initial versions of the modeling code
  - Alternate modeling approaches
  - Additional plotting and visualization tools

  ⚠️ **Files within** ''' bash Master_Data_Sets_FULL''' **may not be stable**

---

## 📧 Points of Contact

- **Corresponding Author:** Ahmed Khan  ahmedk2@illinois.edu
- **Institution:** Department of Aerospace Engineering, University of Illinois at Urbana-Champaign
- **Collaborators:** Tiger Hou  linyhi2@illinois.edu, Seigfried Eggl  eggl@illinois.edu

---

## 💻 Usage

- ✅ **Actively Maintained** – The package is current and stable for general use.
- 🚫  **Not a Standalone Package** – To run the scripts, clone the repository and run them directly with Python. We reccommend using an Anaconda enviroment. 
- ⚠️ **Requires Manual Data Handling** – Inputs will need to be directly inputed by the user.

---

## 🚀 How to Use

### 🧰 Prerequisites

Make sure the following are installed:

```bash
pip install numpy pandas matplotlib scipy lightkurve scipy astropy unpopular scienceplots astroquery 
```

### ▶️ Running the Code

From the root directory:

```bash
python final_code_base.py
```

### 📈 Sample Usage

Here is an example of using the code inside `final_code_base.py` to model a star:

```python
from final_code_base import fit_light_curve, load_data

# Load light curve data for a specific Delta Scuti star
time, flux = load_data("Master_Data_Sets_FULL/KIC_9851822.csv")

# Fit the model
fit_result = fit_light_curve(time, flux, num_terms=5)

# Plot the fit
fit_result.plot()
```

---

## 📌 Documentation & Links

- [Delta Scuti Variable Stars Overview – AAVSO](https://www.aavso.org/delta-scuti-variables)
- [Fourier Analysis in Variable Star Research – NASA ADS](https://ui.adsabs.harvard.edu/)
- (Add links to internal documentation or GitHub repo if available)

---

## 📂 Summary of Files: 

| Files                      | Purpose                        |
| -------------------------- | ------------------------------ |
| `final_code_base.py`       | Core modeling Framework        |
| `Master_Data_Sets_FULL/`   | Complete Datasets              |
| `Legacy_supporting_files/` | Non-essential, legacy scripts  |



