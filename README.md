# 📦 Delta Scuti Light Curve Modeling Framework

This project provides a computational framework for modeling the light curves of **Delta Scuti variable stars**. These stars exhibit periodic brightness variations due to radial and non-radial pulsations, and our tools are designed to analyze and fit these variations based on observational data.

---

## 📚 Package Contents

- ``\
  This is the **core script** that contains all the **functions and classes** used to model, fit, and analyze Delta Scuti light curves. It includes:

  - Fourier decomposition modeling
  - R² and residual analysis
  - Preprocessing and smoothing routines
  - Plotting utilities for light curves

- ``\
  Directory containing the **entire set of analyzed observational data**, including:

  - Cleaned light curve datasets
  - Star metadata (KIC, TIC IDs, pulsation periods, etc.)
  - Sample outputs (e.g., fitted parameters and residuals)

- ``\
  Contains **archived and non-essential scripts**, such as:

  - Initial versions of the modeling code
  - Alternate modeling approaches (e.g., spline fits, Lomb-Scargle periodograms)
  - Visualization tools used during early development

---

## 📞 Points of Contact

- **Lead Developer/PI:** [Your Name or Email]
- **Institution/Organization:** [e.g., University Name / Department]
- **Collaborators:** [List collaborators if applicable]

---

## 📌 Status

- ✅ **Actively Maintained** – The package is current and stable for general use.
- ⚠️ **Not on PyPI / Not a Standalone Package** – To run the scripts, clone the repository and run them directly with Python.
- 📁 **Requires Manual Data Handling** – Input datasets are located in the `Master_Data_Sets_FULL/` directory and must be referenced in your local runtime.

---

## 🚀 How to Use

### 🧰 Prerequisites

Make sure the following are installed:

```bash
pip install numpy pandas matplotlib scipy
```

### ▶️ Running the Code

From the root directory:

```bash
python final_code_base.py
```

Or if integrated into a Bazel workflow:

```bash
bazel run :final_code_base
```

To run tests (if available):

```bash
bazel test :test_final_code_base
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

## 📂 File Summary Table

| File / Folder              | Purpose                        |
| -------------------------- | ------------------------------ |
| `final_code_base.py`       | Core modeling logic            |
| `Master_Data_Sets_FULL/`   | Complete analyzed data archive |
| `Legacy_supporting_files/` | Non-essential, legacy scripts  |



