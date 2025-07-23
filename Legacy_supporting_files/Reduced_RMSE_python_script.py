import numpy as np 
import lightkurve as lk 
import matplotlib as plt
import pandas as pd 
import csv 


def compute(name_of_CSV): 
    try:
            df = pd.read_csv(name_of_CSV)
            
                        
            KIC_list = df['KIC'].dropna().astype(str).tolist()
            MSE_list = df['MSE'].dropna().astype(float).tolist()

    except Exception as e:
        print(f"Error loading TIC IDs: {e}")
        return [] 
    amplitude_vectors = []
    for i in range(len(KIC_list)):
        lc = lk.search_lightcurve(f"KIC {KIC_list[i]}").download_all().stitch().remove_outliers(sigma = 5.0)
        flux =  lc.flux.value
        maxima_amplitude = np.percentile(flux, 95) - 1
        maxima_amplitude *= maxima_amplitude
        print(maxima_amplitude)
        amplitude_vectors.append(maxima_amplitude)
    results = []
    count = 0 
    while(count < len(amplitude_vectors)):
        MSE_list[count] = MSE_list[count]/amplitude_vectors[count]
        MSE_list[count] = np.sqrt(MSE_list[count])
        results.append({'KIC': KIC_list[count], 'RMSE': MSE_list[count]})
        count += 1
         
    df = pd.DataFrame(results)
    df.to_csv('KeplerStarsOutput_corrected_RMSE.csv', index=False)



compute(r'ResearchPython\KeplerStarsOutput_combined - Copy.csv')




