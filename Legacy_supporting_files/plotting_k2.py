import numpy as np
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import pandas as pd
import scienceplots
def plot_aitoff():
    plt.style.use(['science', 'no-latex'])
    plt.rcParams.update({'figure.dpi': '250'})
  
    stars = [
        "EPIC 211945791",
        "EPIC 211115721",
        "EPIC 211044267",
        "EPIC 211914004",
        "Fg Virginis",
        "V1228 Tau",
        "V534 Tau",
        "V624 Tau",
        "V647 Tau",
        "V650 Tau"
    ]

    Simbad.add_votable_fields('ra', 'dec')

    ra_list, dec_list, names = [], [], []

    for star in stars:
        try:
            result = Simbad.query_object(star)
            if result is not None:
                ra = result['ra'][0]
                dec = result['dec'][0]
                ra_list.append(ra)
                dec_list.append(dec)
                names.append(star)
                print(f"Got {star}: RA={ra}, DEC={dec}")
            else:
                print(f"Star {star} not found in SIMBAD.")
        except Exception as e:
            print(f"Error fetching {star}: {e}")
    df = pd.DataFrame({'Name': names,'Right Ascension': ra_list, 'Declination': dec_list})
    df.to_csv(r"C:\Users\ahmed\Downloads\eps_k2\location_ra_dec_new.csv")
    ra_shifted = np.where(np.array(ra_list) > 180, np.array(ra_list) - 360, np.array(ra_list))
    ra_rad = np.radians(ra_shifted)
    dec_rad = np.radians(dec_list)
    plt.figure(figsize=(12, 7))
    ax = plt.subplot(111, projection="aitoff")
    ax.scatter(ra_rad, dec_rad, s=3, color="red", marker="o", label = "K2 Stars")
    
    df = pd.read_csv(r"C:\Users\ahmed\research_delta\ResearchPython\Master_Data_Sets_FULL\KEPLER\KeplerStarsOutput_fixed.csv")
    kep_RA = df['Right Ascension']
    kep_DEC =df['Declination']
    ra_list = kep_RA
    dec_list = kep_DEC
    ra_shifted = np.where(np.array(ra_list) > 180, np.array(ra_list) - 360, np.array(ra_list))
    ra_rad = np.radians(ra_shifted)
    dec_rad = np.radians(dec_list)
    ax.scatter(ra_rad, dec_rad, s=3, marker="o", label = "Kepler Stars")
    xticks_deg = np.arange(-180, 181, 60)
    yticks_deg = np.arange(-90, 91, 30)
    xtick_labels = [f"{int(t)}" if abs(t) != 180 else '' for t in xticks_deg]
    ytick_labels = [f"{int(t)}" if abs(t) != 90 else '' for t in yticks_deg]
    ax.set_xticks(np.radians(xticks_deg))
    ax.set_xticklabels(xtick_labels)
    ax.set_yticks(np.radians(yticks_deg))
    ax.set_yticklabels(ytick_labels)
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.legend(loc='upper left', bbox_to_anchor=(0.75, 1.2)) 

    ax.grid(True)
    plt.show()

plot_aitoff()
