import numpy as np
import matplotlib.pyplot as pt 
import lightkurve as lk 
from astropy.timeseries import LombScargle
from lightkurve import search_lightcurve
from lightkurve import search_targetpixelfile
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.table import Table



vizier = Vizier(columns=["*", "_RAJ2000", "_DEJ2000"])
vizier.ROW_LIMIT = 250
hip_catalog = "I/239/hip_va_1"
hip_result = vizier.query_constraints(catalog=hip_catalog, VarType="DSCT")
hip_table = hip_result[0] if hip_result else Table()
hip_DCEPS_result = vizier.query_constraints(catalog = hip_catalog, VarType = "DCEPS")
hip_DCEPS_table = hip_result[0] if hip_DCEPS_result else Table()
hip_DCEP_result = vizier.query_constraints(catalog = hip_catalog, VarType = "DCEP")
hip_DCEP_table = hip_result[0] if hip_DCEP_result else Table()

print(f"HIPPARCOS Delta Scuti stars (DSCT): {len(hip_table)} entries")
print(hip_result[0])
print("\n")
print(f"HIPPARCOS Delta Scuti stars (DCEPS): {len(hip_DCEPS_table)} entries")
print(hip_DCEPS_result[0])
print("\n")
print(f"HIPPARCOS Delta Scuti stars (DCEP): {len(hip_DCEP_table)} entries")
print(hip_DCEP_result[0])

