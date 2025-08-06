# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:04:19 2024

@author: ahmed
"""

import numpy as np
import matplotlib.pyplot as pt 
import lightkurve as lk 
from astropy.timeseries import LombScargle
from lightkurve import search_lightcurve
from lightkurve import search_targetpixelfile
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, lfilter
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from astroquery.mast import Tesscut
from astroquery.mast.utils import parse_input_location
import matplotlib.pyplot as plt
import unpopular

import scienceplots

def getLightCurveData(nameOfStar):
    search_result = lk.search_lightcurve(nameOfStar, quarter=(6,7,8))
    lc = search_result.download_all()
    time = lc.time.value
    flux = lc.time.value
    pt.plot(time, flux)
    
def getPeriodogramData(nameOfStar): 
    x = lk.search_targetpixelfile(nameOfStar).download()
    y = x.to_lightcurve()
    z = y.to_periodogram()
    z.smooth(method='logmedian', filter_width=0.1).plot(linewidth=2,  color='red', label='Smoothed', scale='log')
    z.plot(scale = 'log')
    #return z

def compGetPeriodogramData(nameOfStar): 
    #x = lk.search_targetpixelfile(nameOfStar).download().to_lightcurve()
    #y = lk.search_lightcurve(nameOfStar, quarter = (6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    
    
    master_flux = []
    master_time = []
    #    i+= 1
    #    continue
    search_result = lk.search_tesscut(target=f"{nameOfStar}")
    search_result = search_result[0]
    print(search_result)
    tpf_collection = search_result.download_all(cutout_size=50)
    for l in tpf_collection:
        s = unpopular.Source(l.path, remove_bad=True)
        s.set_aperture(rowlims=[25, 26], collims=[25, 26])
        s.add_cpm_model(exclusion_size=5, n=64, predictor_method="similar_brightness")
        s.set_regs([0.1])
        s.holdout_fit_predict(k=100);

        aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux")
        aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction", weighting=None)
        #pt.plot(s.time, aperture_normalized_flux, ".", c="k", ms=8, label="Normalized Flux")
        #pt.plot(s.time, aperture_cpm_prediction, "-", lw=3, c="C3", alpha=0.8, label="CPM Prediction")
        #pt.xlabel("Time - 2457000 [Days]", fontsize=30)
        #pt.ylabel("Normalized Flux", fontsize=30)
        #pt.tick_params(labelsize=20)
        #pt.legend(fontsize=30)
        
        apt_detrended_flux = s.get_aperture_lc(data_type="cpm_subtracted_flux")
        min_val = np.percentile(apt_detrended_flux, 5)
        max_val = np.percentile(apt_detrended_flux,95)

        normalized_data = (2 * (apt_detrended_flux - min_val) / (max_val - min_val) - 1) 
        master_flux.extend(normalized_data)
        master_time.extend(s.time)
    
    y = lk.TessLightCurve(time=s.time, flux=apt_detrended_flux) #lk.search_lightcurve(nameOfStar).download().remove_outliers(sigma = 5.0)
    #print(y.author)
    ## quarter term ^^^  y = lk.search_lightcurve(nameOfStar, quarter=(6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    z = y.to_periodogram()
    plt.plot(z.frequency, z.power)
    plt.show()
    
    #z.smooth(method='logmedian', filter_width=0.1).plot(linewidth=2,  color='red', label='Smoothed', scale='log')
    return z, y

def GetProperites(periodogram):
    periodogram.show_properties()
    np.set_printoptions(threshold=np.inf)
    b =  periodogram.frequency
    print(b)
  #  print(np.std(b))

def sine_model(t, amplitude, phase, frequency, offset):
    return amplitude * np.sin(2 * np.pi * frequency * t + phase) + offset

def identifyPeaks(nameOfStar, lowerscalar = 0.2):
    pg, lightc = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * lowerscalar, max_power * 1.1])
    #pt.figure(figsize=(10, 6))
    #pt.plot(pg.frequency, pg.power, label='Periodogram')
    x = pg.frequency[peaks]
    #for i in x:
     #   print(i.value)
    y = pg.power[peaks]
    filtered_peaks = []
    
    for i in range(len(x)):
        if len(filtered_peaks) == 0:
            if(x[i].value >= 1):
                filtered_peaks.append(peaks[i])
        else:
            
            if np.abs(x[i].value - pg.frequency[filtered_peaks[-1]].value) <= 0.3:
                if(y[i].value > pg.power[filtered_peaks[-1]].value):
                    filtered_peaks[-1]= peaks[i]
            else: 
                filtered_peaks.append(peaks[i])
    if(len(pg.frequency[filtered_peaks]) > 10):
        #print("hellll")
        return -1, 0        
    pt.scatter(pg.frequency[filtered_peaks], pg.power[filtered_peaks], color='red', zorder=5, label='Local Maxima')
    pt.xlabel('Frequency (cycles/BKJD)')
    pt.ylabel('Power')
    pt.title('Periodogram with Local Maxima: '+ nameOfStar)
    pt.legend()
    pt.show()
    
    return(pg.frequency[filtered_peaks], lightc, pg.power[filtered_peaks])

def identifyPeaksPowerComp(nameOfStar):
    pg, ltcurves = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * 0.2, max_power * 1.1])
    #pt.figure(figsize=(10, 6))
    #pt.plot(pg.frequency, pg.power, label='Periodogram')
    x = pg.frequency[peaks]
    #for i in x:
     #   print(i.value)
    y = pg.power[peaks]
    filtered_peaks = []
   
    for i in range(len(x)):
        if len(filtered_peaks) == 0:
            if(x[i].value >= 1):
                filtered_peaks.append(peaks[i])
        else:
            
            if np.abs(x[i].value - pg.frequency[filtered_peaks[-1]].value) <= 0.3:
                if(y[i].value > pg.power[filtered_peaks[-1]].value):
                    filtered_peaks[-1]= peaks[i]
            else: 
                filtered_peaks.append(peaks[i])
    if(len(pg.frequency[filtered_peaks]) > 10):
        return -1, 0
    #pt.scatter(pg.frequency[filtered_peaks], pg.power[filtered_peaks], color='red', zorder=5, label='Local Maxima')
    #pt.xlabel('Frequency (cycles/BKJD)')
    #pt.ylabel('Power')
    #pt.title('Periodogram with Local Maxima: '+ nameOfStar)
    #pt.legend()
    #pt.show()
    print(pg.frequency[filtered_peaks])
    print(pg.power[filtered_peaks])
    return(pg.power[filtered_peaks], ltcurves)


def guessHelper(a,bounds1,search_result, frequencyfitted):
    #frequencyfitted, lc = identifyPeaks(a)
    lc = search_result
    #lc.plot()
    #pt.show()
    b = 0 
    c = []
    while b < len(frequencyfitted):
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
        time = lc.time.value
        flux =  lc.flux.value
        vi = np.isfinite(time) & np.isfinite(flux)
        time = time[vi]
        flux = flux[vi]
      
        flux_range = (np.percentile(flux, 95) - np.percentile(flux, 5))
        amplitude_guess = flux_range * 0.5
        phase_guess = 0  
        frequency_guess = frequencyfitted[b].value
        offset_guess = np.mean(flux)
        
        ig = [amplitude_guess, phase_guess, frequency_guess, offset_guess]
        
        # Adding bounds: to force some values
        #bounds = ([0.55*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        bounds = ([0.95 * amplitude_guess, 0, 0.9*frequency_guess, np.min(flux)], [1.05 * amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])

        if len(time) == 0 or len(flux) == 0:
            raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=999999999, method='dogbox')
        amplitude, phase, frequency, offset = params
        
        fit_c = sine_model(time, *params)  
        c.append(fit_c)
        b = b + 1
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return c

def guessLegacy(a,bounds1):
    frequencyfitted = identifyPeaks(a)
    search_result = lk.search_lightcurve(a, quarter=(6,7,8))
    lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
    #lc.plot()
    #pt.show()
    b = 0 
    c = []
    while b < len(frequencyfitted):
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
        time = lc.time.value
        flux =  lc.flux.value
        vi = np.isfinite(time) & np.isfinite(flux)
        time = time[vi]
        flux = flux[vi]
      
        flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
        amplitude_guess = flux_range
        phase_guess = 0  
        frequency_guess = frequencyfitted[b].value
        offset_guess = np.mean(flux)
        
        ig = [amplitude_guess, phase_guess, frequency_guess, offset_guess]
        
        # Adding bounds: to force some values
        bounds = ([0.65*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        #bounds = ([0, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])

        if len(time) == 0 or len(flux) == 0:
            raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=999999999, method='trf')
        amplitude, phase, frequency, offset = params
        
        fit_c = sine_model(time, *params)  
        #pt.plot(time, fit_c, 'o-', )
        c.append(fit_c)
        b = b + 1
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return c

def guessIterative(a,bound):
    frequencyfitted = identifyPeaks(a)
    search_result = lk.search_lightcurve(a, quarter=(6,7,8))
    lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
    lc.plot()
    pt.show()
    b = 0 
    c = []
    while b < len(frequencyfitted):
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
        time = lc.time.value
        flux =  lc.flux.value
        vi = np.isfinite(time) & np.isfinite(flux)
        time = time[vi]
        flux = flux[vi]
      
        flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
        amplitude_guess = flux_range
        phase_guess = 0  
        frequency_guess = frequencyfitted[b].value
        offset_guess = np.mean(flux)
        
        ig = [amplitude_guess, phase_guess, frequency_guess, offset_guess]
        
        # Adding bounds: to force some values
        #bounds = ([0.55*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        bounds = ([0.55**amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])

        if len(time) == 0 or len(flux) == 0:
            raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=999999999, method='trf')
        amplitude, phase, frequency, offset = params
        fit_c = sine_model(time, *params)
        #residuals = flux - (fit_c)
        #meanSquare = np.sum((residuals)**2)/len(flux)
        #if(meanSquare < bestfit):
        #    bestfitSineCurve = fit_c
        #    bestfit = meanSquare
        #    print(a)
        c.append(fit_c)
        b = b + 1
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return c

def guessActual(a):
    frequencyfitted, search_result, powers = identifyPeaks(a)
    lc = search_result
    #lc.plot()
    #pt.show()
    b = 0 
    c = []
    params_list = []
    time = lc.time.value
    flux =  lc.flux.value
    print(f"need to iterate: +  {len(frequencyfitted)} + times")
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
    amplitude_guess = flux_range
    phase_guess = 0 
    offset_guess = np.mean(flux)
    while b < len(frequencyfitted):
        
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
 
        frequency_guess = frequencyfitted[b].value
        ig = [0.75*amplitude_guess, phase_guess, frequency_guess, offset_guess]
        # Adding bounds: to force some values of amplitude
        bounds = ([0.55*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
        #bounds = ([y*ampliude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        

        if len(time) == 0 or len(flux) == 0:
              raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, method='dogbox', maxfev = 9999999)
        amplitude, phase, frequency, offset = params
        fit_c = sine_model(time, *params)
        #amplitude, phase, frequency, offset = params
        #residuals = flux - (fit_c)
        #meanSquare = np.sum((residuals)**2)/len(flux)
        c.append(fit_c)
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1
        #flux -= fit_c
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return params_list, lc 

def guessActual_refined(a):
    frequencyfitted, search_result, powers = identifyPeaks(a)
    lc = search_result
    #lc.plot()
    #pt.show()
    b = 0 
    c = []
    params_list = []
    time = lc.time.value
    flux =  lc.flux.value
    flux_org = flux
    print(f"need to iterate: +  {len(frequencyfitted)} + times")
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
     
    total_power = np.sum(powers)
    while b < len(frequencyfitted):
        flux = flux_org * (powers[b]/total_power)
        offset_guess = np.mean(flux)
        flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
        amplitude_guess = flux_range
        phase_guess = 0
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
 
        frequency_guess = frequencyfitted[b].value
        amplitude_guess_scale = 0.75
        ig = [amplitude_guess_scale*amplitude_guess, phase_guess, frequency_guess, offset_guess]
        # Adding bounds: to force some values of amplitude
        amplitude_scale = 0.5
        bounds = ([amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
        #bounds = ([y*ampliude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
         

        if len(time) == 0 or len(flux) == 0:
              raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, method='dogbox')
        amplitude, phase, frequency, offset = params
        fit_c = sine_model(time, *params)
        bestmean = getResiduals(fit_c, flux)
        bestFitAchieved = False

        while not bestFitAchieved:
            low_amplitude_scale = amplitude_scale * 0.9
            high_amplitude_scale = amplitude_scale * 1.1
            low_amplitude_guess_scale = low_amplitude_scale * 1.1
            high_amplitude_guess_scale = high_amplitude_scale * 1.1
            
            
            ig_low = [low_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            ig_high = [high_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            
            bounds_low = ([low_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], 
                        [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.percentile(flux,95)])
            bounds_high = ([high_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], 
                        [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.percentile(flux,95)])
            
           
            try:
                params_low, _ = curve_fit(sine_model, time, flux, p0=ig_low, bounds=bounds_low, method='dogbox')
                fit_low = sine_model(time, *params_low)
                fit_low_MSE = getResiduals(fit_low, flux)
            except:
                fit_low_MSE = np.inf
                
            try:
                params_high, _ = curve_fit(sine_model, time, flux, p0=ig_high, bounds=bounds_high, method='dogbox')
                fit_high = sine_model(time, *params_high)
                fit_high_MSE = getResiduals(fit_high, flux)
            except:
                fit_high_MSE = np.inf
            
            current_MSE = bestmean
            
            
            options = [
                ('low', fit_low_MSE, params_low, low_amplitude_scale),
                ('current', current_MSE, params, amplitude_scale),
                ('high', fit_high_MSE, params_high, high_amplitude_scale)
            ]
            
            
            best_option = min(options, key=lambda x: x[1])
            
            if best_option[0] == 'current':
                bestFitAchieved = True
                print("-1")
            else:
                
                bestmean = best_option[1]
                params = best_option[2]
                amplitude_scale = best_option[3]
                fit_c = sine_model(time, *params)
                print(f"Switching to {best_option[0]} fit")
        """
        while(bestFitAchieved == False): 
            low_amplitude_scale = amplitude_scale*  0.9
            high_amplitude_scale =  amplitude_scale*  1.1
            low_amplitude_guess_scale = low_amplitude_scale * 1.1
            high_amplitude_guess_scale = high_amplitude_scale * 1.1
            ig_low = [low_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            ig_high =  [high_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            bounds_high = ([high_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [ amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
            bounds_low = ([low_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [ amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
            
            params_low, _ = curve_fit(sine_model, time, flux, p0=ig_low, bounds=bounds_low, method='dogbox')
            params_high, _ = curve_fit(sine_model, time, flux, p0=ig_high, bounds=bounds_high, method='dogbox')
            fit_low = sine_model(time, *params_low)
            fit_high = sine_model(time, *params_high)
            fit_low_MSE = getResiduals(fit_low, flux)
            fit_high_MSE = getResiduals(fit_high, flux)
            order_fit = np.array([fit_low, fit_high, fit_c])
            order_arg = np.array([fit_low_MSE, fit_high_MSE, bestmean])
            order_params = np.array([params_low, params_high])
            #order_scales_guess = np.array([low_amplitude_guess_scale, high_amplitude_guess_scale])
            order_scales = np.array([low_amplitude_scale, high_amplitude_scale])
            best_index = np.argmin(order_arg)
            minimum_val = np.min(order_arg)
            if best_index != 2:
                print(best_index)
                bestmean = order_arg[best_index]
                amplitude_scale = order_scales[best_index]
                fit_c = order_fit[best_index]
                params = order_params[best_index]
            else:
              bestFitAchieved = True
              print("-1")
              break
            """

            

        #ig_refined = [amplitude + 0.00001, phase, frequency, offset]
        #amplitude, phase, frequency, offset = params
        #residuals = flux - (fit_c)
        #meanSquare = np.sum((residuals)**2)/len(flux)
        #bounds_refined =  ([amplitude, -2*np.pi, frequency * 0.8, np.percentile(flux,5)], [amplitude * 1.5, 2*np.pi, frequency * 1.2,  np.percentile(flux,95)])
        #params_refined, _ = curve_fit(sine_model, time, flux, p0=ig_refined, bounds=bounds_refined, method='dogbox')
        #amplitude, phase, frequency, offset = params
        #fit_refined = sine_model(time, *params_refined)
        c.append(fit_c)

        amplitude, phase, frequency, offset = params
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1
        
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return params_list, lc 
#params_list.append((amplitude, phase, frequency, offset)) 

def guessActual_refined_second_iteration(a, scalar, frequencyfitted, search_result, powers):
    lc = search_result
    #lc.plot()
    #pt.show()
    b = 0 
    c = []
    params_list = []
    time = lc.time.value
    flux =  lc.flux.value
    flux_org = flux
    print(f"need to iterate: +  {len(frequencyfitted)} + times")
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    if(len(frequencyfitted) == 0):
        return [0],[0],[0]

    while b < len(frequencyfitted):
        offset_guess = np.mean(flux)
        flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
        amplitude_guess = flux_range
        amplitude_guess_upper = amplitude_guess
        phase_guess = 0
        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
 
        frequency_guess = frequencyfitted[b].value
        amplitude_guess_scale =  scalar * 1.5
        ig = [amplitude_guess_scale*amplitude_guess, phase_guess, frequency_guess, offset_guess]
        # Adding bounds: to force some values of amplitude
        amplitude_scale = scalar
        if(amplitude_guess_scale>= 1):
            amplitude_guess_upper *= (amplitude_guess_scale + 0.2)

        bounds = ([amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [amplitude_guess_upper, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
        #bounds = ([y*ampliude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        

        if len(time) == 0 or len(flux) == 0:
              raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, method='dogbox', maxfev=999999999)
        amplitude, phase, frequency, offset = params
        fit_c = sine_model(time, *params)
        
        """
        while(bestFitAchieved == False): 
            low_amplitude_scale = amplitude_scale*  0.9
            high_amplitude_scale =  amplitude_scale*  1.1
            low_amplitude_guess_scale = low_amplitude_scale * 1.1
            high_amplitude_guess_scale = high_amplitude_scale * 1.1
            ig_low = [low_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            ig_high =  [high_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            bounds_high = ([high_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [ amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
            bounds_low = ([low_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], [ amplitude_guess, 2*np.pi, 1.1*frequency_guess,  np.percentile(flux,95)])
            
            params_low, _ = curve_fit(sine_model, time, flux, p0=ig_low, bounds=bounds_low, method='dogbox')
            params_high, _ = curve_fit(sine_model, time, flux, p0=ig_high, bounds=bounds_high, method='dogbox')
            fit_low = sine_model(time, *params_low)
            fit_high = sine_model(time, *params_high)
            fit_low_MSE = getResiduals(fit_low, flux)
            fit_high_MSE = getResiduals(fit_high, flux)
            order_fit = np.array([fit_low, fit_high, fit_c])
            order_arg = np.array([fit_low_MSE, fit_high_MSE, bestmean])
            order_params = np.array([params_low, params_high])
            #order_scales_guess = np.array([low_amplitude_guess_scale, high_amplitude_guess_scale])
            order_scales = np.array([low_amplitude_scale, high_amplitude_scale])
            best_index = np.argmin(order_arg)
            minimum_val = np.min(order_arg)
            if best_index != 2:
                print(best_index)
                bestmean = order_arg[best_index]
                amplitude_scale = order_scales[best_index]
                fit_c = order_fit[best_index]
                params = order_params[best_index]
            else:
              bestFitAchieved = True
              print("-1")
              break
            """

            

        #ig_refined = [amplitude + 0.00001, phase, frequency, offset]
        #amplitude, phase, frequency, offset = params
        #residuals = flux - (fit_c)
        #meanSquare = np.sum((residuals)**2)/len(flux)
        #bounds_refined =  ([amplitude, -2*np.pi, frequency * 0.8, np.percentile(flux,5)], [amplitude * 1.5, 2*np.pi, frequency * 1.2,  np.percentile(flux,95)])
        #params_refined, _ = curve_fit(sine_model, time, flux, p0=ig_refined, bounds=bounds_refined, method='dogbox')
        #amplitude, phase, frequency, offset = params
        #fit_refined = sine_model(time, *params_refined)
        c.append(fit_c)

        amplitude, phase, frequency, offset = params
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1
        print("passed")
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return params_list, lc, c
#params_list.append((amplitude, phase, frequency, offset)) 



def align_arrays(time, flux):
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    return time, flux

def getMeanSquaredResidual(a, search_result, frequency, powerofpeaks_arg):
        bestmeanSquare = 100000
        bestBound = 0
        lc = search_result
        for bounds1 in range(54,56): #######
            listofsines = guessHelper(a,bounds1, search_result, frequency)
            addedTogether = 0
            time = lc.time.value
            flux = lc.flux.value
            #flux_err = np.mean(lc.flux_err.value)
            time, flux = align_arrays(time,flux)
            powerOfPeaks = powerofpeaks_arg
            powerOfPeaks = powerOfPeaks.value
            p = 0 
            total_weight = 0 
            while (p < len(powerOfPeaks)):
                sinInterpolated = interpolate(time, listofsines[p], time)
                weight = powerOfPeaks[p]
                #print(weight) 
                total_weight += weight
                addedTogether += (weight * sinInterpolated)
                p += 1
                addedTogether  = addedTogether/total_weight
                residuals = flux - (addedTogether)
                meanSquare = np.sum((residuals)**2)/len(flux)
                if(meanSquare < bestmeanSquare):
                    bestmeanSquare = meanSquare
                    bestBound = bounds1
                    #print(meanSquare)
                    #print(bounds1/100)
        #al = len(flux)  
        #p = 4 * len(listofsines)  
        #reduced_chi_squared = chi_squared / (al - p)
        #return reduced_chi_squared
        print(bestmeanSquare)
        return bestmeanSquare, bestBound/100


def getResiduals(fit, flux): 
    residuals = flux - fit
    meanSquare = np.sum((residuals)**2)/len(flux)
    return meanSquare


def getCompositeSine(a):
        listofsines = guessLegacy(a,0)
        addedTogether = 0
        search_result = lk.search_lightcurve(a,quarter=(6,7,8))
        lc = search_result.download_all().stitch()
        time = lc.time.value
        flux = lc.flux.value
        time, flux = align_arrays(time,flux)
        powerOfPeaks = identifyPeaksPowerComp(a).value


        p = 0 
        total_weight = 0 
        sine_print_terms = []
        while (p < len(powerOfPeaks)):
           #amplitude, phase, frequency, offset = listofsines[p]
           sinInterpolated = interpolate(time, listofsines[p], time)
           weight = powerOfPeaks[p]
           print(weight) 
           total_weight += weight
           addedTogether += (weight * sinInterpolated)
           #sine_print_terms.append(f"{amplitude:.2f} * sin(2π * {frequency:.2f} * t + {phase:.2f})")
           p += 1
        addedTogether  = addedTogether/total_weight
        #print(f"Composite Sine Function for {a}:")
        #print("f(t) = " + " + ".join(sine_print_terms))
        #print(total_weight)
        return addedTogether

def getCompositeSine2(a):
        powerOfPeaks, _ = identifyPeaksPowerComp(a)
        print(len(powerOfPeaks))
        powerOfPeaks = powerOfPeaks.value
        listofsines, lc = guessActual(a)
        addedTogether = 0
        time = lc.time.value
        flux = lc.flux.value
        time, flux = align_arrays(time,flux)
        p = 0 
        total_weight = np.sum(powerOfPeaks)
        sine_print_terms = []
        while (p < len(powerOfPeaks)):
           
            amplitude, phase, frequency, offset = listofsines[p]
            sinInterpolated = amplitude * np.sin(2 * np.pi * frequency * time + phase) + offset
            weight = powerOfPeaks[p]  
            amplitude = amplitude * (weight/total_weight)
            offset = offset * (weight/total_weight)
            addedTogether += (weight/total_weight) * sinInterpolated
            #addedTogether += sinInterpolated
            sine_print_terms.append(f"{amplitude:.4f} * sin(2π * {frequency:.4f} * t + {phase:.4f}) + {offset:.4f}")
            p += 1
        #addedTogether  = addedTogether/total_weight
        print(f"Composite Sine Function for {a}:")
        print("f(t) = " + " + ".join(sine_print_terms))
        print(total_weight)
        return addedTogether, lc     


def getCompositeSine2_second_test(a):
        powerOfPeaks, _ = identifyPeaksPowerComp(a)
        #if(powerOfPeaks == [-1]):
        #    
        #    return [-10], 0,"0"
        print(len(powerOfPeaks))
        powerOfPeaks = powerOfPeaks.value
        frequencyfitted2, search_result2, powers2 = identifyPeaks(a)
        amplitude_scale = 0.5
        listofsines, lc, _ = guessActual_refined_second_iteration(a, amplitude_scale, frequencyfitted2, search_result2, powers2)
        if(listofsines == [0]):
            return [-10],[0],"0"
        listofindexs =[]
        addedTogether = 0
        time = lc.time.value
        flux = lc.flux.value
        time, flux = align_arrays(time,flux)
        p = 0 
        total_weight = np.sum(powerOfPeaks)
        sine_print_terms = []

        
        while (p < len(listofsines)):
           
            amplitude, phase, frequency, offset = listofsines[p]
            sinInterpolated = amplitude * np.sin(2 * np.pi * frequency * time + phase) + offset
            weight = powerOfPeaks[p]  
            amplitude = amplitude * (weight/total_weight)
            offset = offset * (weight/total_weight)
            addedTogether += (weight/total_weight) * sinInterpolated
            p+=1
        bestmean = getResiduals(addedTogether, flux)
        bestFitAchieved = False
        while(bestFitAchieved == False): 
            low_amplitude_scale = amplitude_scale*  0.9
            high_amplitude_scale =  amplitude_scale*  1.1
            
            lower, _,fits_low = guessActual_refined_second_iteration(a, low_amplitude_scale, frequencyfitted2, search_result2, powers2)
            upper, _, fits_high = guessActual_refined_second_iteration(a, high_amplitude_scale, frequencyfitted2, search_result2, powers2)
            lowertot = 0
            uppertot = 0 
            countinner = 0
            while (countinner < len(listofsines)):
           
                amplitude, phase, frequency, offset = lower[countinner]
                sinInterpolated = amplitude * np.sin(2 * np.pi * frequency * time + phase) + offset
                weight = powerOfPeaks[countinner]  
                amplitude = amplitude * (weight/total_weight)
                offset = offset * (weight/total_weight)
                lowertot += (weight/total_weight) * sinInterpolated

                amplitude, phase, frequency, offset = upper[countinner]
                sinInterpolated = amplitude * np.sin(2 * np.pi * frequency * time + phase) + offset
                weight = powerOfPeaks[countinner]  
                amplitude = amplitude * (weight/total_weight)
                offset = offset * (weight/total_weight)
                uppertot += (weight/total_weight) * sinInterpolated

                countinner+=1
            print(lowertot, uppertot)
            fit_low_MSE = getResiduals(lowertot, flux)
            fit_high_MSE = getResiduals(uppertot, flux)
            order_fit = np.array([lowertot, uppertot, addedTogether])
            order_arg = np.array([fit_low_MSE, fit_high_MSE, bestmean])
            order_params = np.array([lower, upper])
            #order_scales_guess = np.array([low_amplitude_guess_scale, high_amplitude_guess_scale])
            order_scales = np.array([low_amplitude_scale, high_amplitude_scale])
            best_index = np.argmin(order_arg)
            minimum_val = np.min(order_arg)
            listofindexs.append(best_index)
            if (best_index < 2):
                print(best_index)
                bestmean = order_arg[best_index]
                amplitude_scale = order_scales[best_index]
                listofsines = order_params[best_index]
                addedTogether = order_fit[best_index]
            else:
              bestFitAchieved = True
              print("-1")
              break
        count = 0 
        newaddedtogether = 0
        while (count < len(listofsines)):
           
            amplitude, phase, frequency, offset = listofsines[count]
            sinInterpolated = amplitude * np.sin(2 * np.pi * frequency * time + phase) + offset
            weight = powerOfPeaks[count]  
            amplitude = amplitude * (weight/total_weight)
            offset = offset * (weight/total_weight)
            newaddedtogether += (weight/total_weight) * sinInterpolated
            #addedTogether += sinInterpolated
            sine_print_terms.append(f"{amplitude:.4f} * sin(2π * {frequency:.4f} * t + {phase:.4f}) + {offset:.4f}")
            count += 1
        #addedTogether  = addedTogether/total_weight
        composite_string  = "f(t) = " + " + ".join(sine_print_terms)
        print(f"Composite Sine Function for {a}:")
        print(composite_string)
        print(total_weight)
        print(listofindexs)
        return newaddedtogether, lc, composite_string    

def plotsidebysideactual_manual(a):
    function, lc = getCompositeSine2(a)
    flux = lc.flux.value
    print(flux)
    print(function)
    time = lc.time.value
    min_length = min(len(flux), len(function))
    flux = flux[:min_length]
    time = time[:min_length]
    function = function[:min_length]
    residuals = flux - function
    print(f"MSE: {np.sum((residuals)**2)/len(flux)}")
    #a = 0.00219*np.sin(2*np.pi*10.33759*time+-0.21704)+ 0.54456 + 0.00183*np.sin(2*np.pi*12.47142*time+-6.28319) + 0.45546
    print(residuals)
    pt.plot(time, residuals, 'o-', color='blue', label='O-C (Observed - Calculated)')
    pt.plot(time, flux, 'o-', color='red', label='Light Curve')
    pt.plot(time, function, 'o-', color='green', label='Curve Fit')
    #pt.plot(time, a, 'o-', color = 'blue')
    pt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    pt.title("O-C Diagram " + str(a))
    pt.xlabel("Time (Days)")
    pt.ylabel("O-C (Flux Difference)")
    pt.legend()
    pt.grid()
    pt.tight_layout()
    pt.show()
    #pt.plot(flux, 'b')
    #pt.plot(function, 'r')
    #pt.plot(residuals, 'g')
    #pt.show()




def plotsidebysideactual(a):
    function, lc, _ = getCompositeSine2_second_test(a)
    
    flux = lc.flux.value
    #print(flux)
    #print(function)
    time = lc.time.value
    min_length = min(len(flux), len(function))
    flux = flux[:min_length]
    time = time[:min_length]
    function = function[:min_length]
    residuals = flux - function
    print(f"MSE: {np.sum((residuals)**2)/len(flux)}")
    #a = 0.00219*np.sin(2*np.pi*10.33759*time+-0.21704)+ 0.54456 + 0.00183*np.sin(2*np.pi*12.47142*time+-6.28319) + 0.45546
    print(residuals)
    pt.plot(time, residuals, 'o-', color='blue', label=' Residuals (Observed - Calculated)')
    pt.plot(time, flux, 'o-', color='red', label='Light Curve')
    pt.plot(time, function, 'o-', color='green', label='Curve Fit')
    #pt.plot(time, a, 'o-', color = 'blue')
    pt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    pt.title("Diagram For " + str(a))
    pt.xlabel("Time -2454833 [BKJD Days]")
    pt.ylabel("Normalized Flux")
    pt.legend()
    pt.grid()
    pt.tight_layout()
    pt.show()
    #pt.plot(flux, 'b')
    #pt.plot(function, 'r')
    #pt.plot(residuals, 'g')
    #pt.show()

def plotsidebyside2(a):
    function = getCompositeSine(a)
    lc = lk.search_lightcurve(a,quarter=(6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    flux = lc.flux.value
    print(flux)
    print(function)
    time = lc.time.value
    min_length = min(len(flux), len(function))
    flux = flux[:min_length]
    time = time[:min_length]
    function = function[:min_length]
    residuals = flux - function
    print(np.sum((residuals)**2)/len(flux))
    #a = 0.00219*np.sin(2*np.pi*10.33759*time+-0.21704)+ 0.54456 + 0.00183*np.sin(2*np.pi*12.47142*time+-6.28319) + 0.45546
    print(residuals)
    pt.plot(time, residuals, 'o-', color='blue', label='O-C (Observed - Calculated)')
    pt.plot(time, flux, 'o-', color='red', label='Light Curve')
    pt.plot(time, function, 'o-', color='green', label='Curve Fit')
    #pt.plot(time, a, 'o-', color = 'blue')
    pt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    pt.title("O-C Diagram " + str(a))
    pt.xlabel("Time (Days)")
    pt.ylabel("O-C (Flux Difference)")
    pt.legend()
    pt.grid()
    pt.tight_layout()
    pt.show()
    #pt.plot(flux, 'b')
    #pt.plot(function, 'r')
    #pt.plot(residuals, 'g')
    #pt.show()

def interpolate(time, flux, target_time):
    ip = interp1d(time, flux, kind='nearest', bounds_error=False, fill_value='extrapolate')
    interpolated_flux = ip(target_time)
    return interpolated_flux 

def guessIteration(a):
    x = 0 
    search_result = lk.search_lightcurve(a, quarter=(6,7,8))
    c = []
    while x < len(search_result):
       currentchi = guessHelper(a, x)
       c.append(currentchi)
       x = x + 1
    #print(f"Reduced Chi-squared Mean: {np.mean(c):.3f}")
    return (np.mean(c))

def seriesOfStars(z): 
    ab = []
    for x in z:
        a = getMeanSquaredResidual(x)
        #plotsidebyside(x)
        print(str(x) + " " + str(a))
        #plotsidebyside(x)
        ab.append(a)
    print("The Best Star is " + str(z[np.argmin(ab)]) + " " + str(np.min(ab)))

def getInfomation(listofStars):
    for x  in listofStars:
        b = lk.search_lightcurve(x, quarter=(6,7,8)).target_name
        print(b.target_name)

def identifyPeaksOfLightcurves(nameofStar,startingTime): 
    Composite_function, lc, _ = getCompositeSine2_second_test(nameofStar) #2 # actual 
    time = lc.time.value
    print(time)
    flux = lc.flux.value
    time, flux = align_arrays(time,flux)
    gradientOfComposite = np.gradient(Composite_function,time)
    gradientOfLightcurve = np.gradient(flux, time)
    print(gradientOfComposite)
    print(gradientOfLightcurve)
    def gettimestamps(gradientarray,time):
        a = []
        b = []
        x = 0 
        lasttimestamp = 0 
        while x < len(gradientarray):
            if (np.absolute(gradientarray[x]) < 0.1): 
                if (lasttimestamp != time[x-1]):
                    a.append(time[x])
                    b.append(flux[x])
                    lasttimestamp = time[x]
                    x +=1
            x += 1
        return a,b
    peaksofcomposite, fluxvalues = gettimestamps(gradientOfComposite, time)
    peaksoflightcurve,fluxvaluesLight = gettimestamps(gradientOfLightcurve, time)
    min_length = min(len(peaksofcomposite), len(peaksoflightcurve))
    peaksofcomposite = np.array(peaksofcomposite[:min_length])
    peaksoflightcurve = np.array(peaksoflightcurve[:min_length])
    print(peaksofcomposite)
    #peaksofcomposite = np.array(peaksofcomposite)[composite_time_mask]
    #peaksoflightcurve = np.array(peaksoflightcurve)[composite_time_mask]
    #pt.scatter(peaksofcomposite, np.array(fluxvalues[:min_length]), color='blue', label='composite')
    distances = cdist(peaksofcomposite.reshape(-1, 1), peaksoflightcurve.reshape(-1, 1))
    indices = np.argmin(distances, axis=1)
    matched_lightcurve_peaks = peaksoflightcurve[indices]
    #np.set_printoptions(threshold=np.inf)
    #print(peaksofcomposite)
    #print(matched_lightcurve_peaks)
    matched_lightcurve_peaks_snipped = matched_lightcurve_peaks
    peaksofcomposite_snipped = peaksofcomposite
    print(f"Snipped peaks lengths: {len(peaksofcomposite_snipped)}, {len(matched_lightcurve_peaks_snipped)}")
    pt.figure(figsize=(12, 8))
    #pt.scatter(peaksofcomposite_snipped, np.ones(peaksofcomposite_snipped.shape[0]), color='blue', label='composite')
    residuals = np.abs(peaksofcomposite_snipped - matched_lightcurve_peaks_snipped)
    #print(residuals.shape)
    #pt.scatter(matched_lightcurve_peaks, np.array(fluxvaluesLight[:min_length]), color='red', label='Light Curve')
    #pt.scatter(matched_lightcurve_peaks_snipped, np.ones(matched_lightcurve_peaks_snipped.shape[0]), color='red', label='Light Curve')
    #pt.plot(matched_lightcurve_peaks_snipped,residuals,"o-", color = "black")
    #print(residuals)
    print(f" The average residual in days is {np.mean(residuals)}")
    #pt.plot(time, flux[:min_length2], 'o-', color='black', label='Light Curve')
    #pt.plot(time, Composite_function[:min_length2], 'o-', color='green', label='Curve Fit')
    pt.tight_layout()
    #pt.draw()
    #pt.show()
    #print("Original peaksofcomposite:", peaksofcomposite[:10])  # first 10 values
    #print("Snipped peaksofcomposite:", peaksofcomposite_snipped[:10])  # first 10 values after snipping
    #print("Original matched_lightcurve_peaks:", peaksoflightcurve[:10])
    #print("Snipped matched_lightcurve_peaks:", matched_lightcurve_peaks_snipped[:10])
    #print("StartingTime (index to cut):", startingTime)
    return np.mean(residuals)


def identifyPeaksOfLightcurves_manual(nameofStar,startingTime): 
    Composite_function, lc = getCompositeSine2(nameofStar) #2 # actual 
    time = lc.time.value
    print(time)
    flux = lc.flux.value
    time, flux = align_arrays(time,flux)
    gradientOfComposite = np.gradient(Composite_function,time)
    gradientOfLightcurve = np.gradient(flux, time)
    print(gradientOfComposite)
    print(gradientOfLightcurve)
    def gettimestamps(gradientarray,time):
        a = []
        b = []
        x = 0 
        lasttimestamp = 0 
        while x < len(gradientarray):
            if (np.absolute(gradientarray[x]) < 0.1): 
                if (lasttimestamp != time[x-1]):
                    a.append(time[x])
                    b.append(flux[x])
                    lasttimestamp = time[x]
                    x +=1
            x += 1
        return a,b
    peaksofcomposite, fluxvalues = gettimestamps(gradientOfComposite, time)
    peaksoflightcurve,fluxvaluesLight = gettimestamps(gradientOfLightcurve, time)
    residuals_flux = flux - Composite_function
    print(f"MSE: {np.sum((residuals_flux)**2)/len(flux)}")
    min_length = min(len(peaksofcomposite), len(peaksoflightcurve))
    peaksofcomposite = np.array(peaksofcomposite[:min_length])
    peaksoflightcurve = np.array(peaksoflightcurve[:min_length])
    print(peaksofcomposite)
    #peaksofcomposite = np.array(peaksofcomposite)[composite_time_mask]
    #peaksoflightcurve = np.array(peaksoflightcurve)[composite_time_mask]
    #pt.scatter(peaksofcomposite, np.array(fluxvalues[:min_length]), color='blue', label='composite')
    distances = cdist(peaksofcomposite.reshape(-1, 1), peaksoflightcurve.reshape(-1, 1))
    indices = np.argmin(distances, axis=1)
    matched_lightcurve_peaks = peaksoflightcurve[indices]
    #np.set_printoptions(threshold=np.inf)
    #print(peaksofcomposite)
    #print(matched_lightcurve_peaks)
    matched_lightcurve_peaks_snipped = matched_lightcurve_peaks
    peaksofcomposite_snipped = peaksofcomposite
    print(f"Snipped peaks lengths: {len(peaksofcomposite_snipped)}, {len(matched_lightcurve_peaks_snipped)}")
    pt.figure(figsize=(12, 8))
    pt.scatter(peaksofcomposite_snipped, np.ones(peaksofcomposite_snipped.shape[0]), color='blue', label='Composite function')
    residuals = np.abs(peaksofcomposite_snipped - matched_lightcurve_peaks_snipped)
    #print(residuals.shape)
    #pt.scatter(matched_lightcurve_peaks, np.array(fluxvaluesLight[:min_length]), color='red', label='Light Curve')
    pt.scatter(matched_lightcurve_peaks_snipped, np.ones(matched_lightcurve_peaks_snipped.shape[0]), color='red', label='Light curve')
    pt.plot(matched_lightcurve_peaks_snipped,residuals,"o-", color = "black", label = "ϵ value")
    #print(residuals)
    print(f" The average residual in days is {np.mean(residuals)}")
    #pt.plot(time, flux[:min_length2], 'o-', color='black', label='Light Curve')
    #pt.plot(time, Composite_function[:min_length2], 'o-', color='green', label='Curve Fit')
    pt.tight_layout()
    pt.legend()
    pt.title("Peak and Trough Time Difference for KIC 3123138")
    pt.xlabel("Time (BKJD days)")
    pt.ylabel("ϵ value (BKJD days)")
    #pt.draw()
    pt.show()
    #print("Original peaksofcomposite:", peaksofcomposite[:10])  # first 10 values
    #print("Snipped peaksofcomposite:", peaksofcomposite_snipped[:10])  # first 10 values after snipping
    #print("Original matched_lightcurve_peaks:", peaksoflightcurve[:10])
    #print("Snipped matched_lightcurve_peaks:", matched_lightcurve_peaks_snipped[:10])
    #print("StartingTime (index to cut):", startingTime)
    return np.mean(residuals)


def get_epsilon_value(star_name, sine_string):

    search_result = lk.search_lightcurve(f"KIC {star_name}")
    lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
    t = lc.time.value
    #sine_string = "0.0020 * np.sin(2* np.pi * (10.3376) * t + -0.2050) + 1 + 0.0017 * np.sin(np.pi * 2 * (12.4714) * t + -6.2832)"
    #sine_string = "0.0020 sin(2π(10.3376)t + -0.2050) + 1 + 0.0017 sin(2π(12.4714)t + -6.2832)"
    #0.0020  * np.sin(2 * np.pi * (10.3376) * t  + -0.2050) + 1 + 0.0017  * np.sin(2 * np.pi * (12.4714) * t  + -6.2832)
    sine_string = sine_string.replace('sin', 'np.sin')
    sine_string = sine_string.replace('2π', '2 * np.pi ')
    #sine_string = sine_string.replace('t', ' * t ')
    sine_string = sine_string.replace("f(t) = ", "")
    #model_profile = eval(sine_string)
    OFFSET = 0
    expected_cadence = 1800  # seconds

    
    def create_model_function(sine_string):
        """Create a callable function from the sine string"""
        def model(t, dt, *params):

            shifted_t = t + dt #+ (OFFSET-2457000)
            print(sine_string)
            return eval(sine_string.replace('t', 'shifted_t'))
        return model

    
    #sine_string = "0.0020 * np.sin(2* np.pi * (10.3376) * t + -0.2050) + 1 + 0.0017 * np.sin(np.pi * 2 * (12.4714) * t + -6.2832)"
    profile_func = create_model_function(sine_string)

    mask = (np.isfinite(lc.flux.value.unmasked))
    all_flux = lc.flux.value[mask]
    all_time = lc.time.value[mask]

    true_time = []
    est_time = []
    t_step = 0.1  
    dt = 1.0  
    t_prev = -np.inf

    for t in all_time:
        if t - t_prev > t_step:
            t_prev = t
        else:
            continue

        
        start_bxjd = t_prev
        mask = (all_time >= start_bxjd) & (all_time <= (start_bxjd + dt))
        time = all_time[mask]
        flux = all_flux[mask]

       
        if len(time) == 0:
            continue
        if abs(time[-1] - start_bxjd - dt) > expected_cadence/86400:
            continue
        if abs(time[0] - start_bxjd) > expected_cadence/86400:
            continue
        if np.any(np.diff(time) > expected_cadence/86400):
            continue

        t_zeroed = time - time[0]
        
        
        t0_list = np.arange(-4,4) * 0.1 + time[0] + np.random.normal(0.01, 0.05)
        t_est_list = []
        
        for t0 in t0_list:
            try:
                popt, pcov = curve_fit(
                    profile_func,
                    xdata=t_zeroed,
                    ydata=flux,
                    p0=t0,
                    xtol=1e-12,
                    maxfev=1000
                )
                t_est_list.append(popt[0])
            except RuntimeError:
                continue
        
        if len(t_est_list) > 0:
            t_est = t_est_list[np.argmin(np.abs(t_est_list-time[0]))]
            true_time.append(time[0])
            est_time.append(t_est)

   
    true_time = np.array(true_time)
    est_time = np.array(est_time)

    
    time_diff = np.diff(true_time)
    mask = time_diff > 3000
    gap_indices = np.where(mask)[0]
    segments = np.split(true_time, gap_indices+1)
    """
    # plotting
    tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5)/100)*100)
    margin = 0.5  # days
    
    fig_oc, axs = pt.subplots(1, len(segments), figsize=(8,3), sharey=True, 
                            gridspec_kw={'wspace': 0, 'hspace': 0},
                            width_ratios=[seg[-1]-seg[0] + margin*2 for seg in segments])

    if not isinstance(axs, np.ndarray):
        axs = [axs]

    for ii, ax in enumerate(axs):
        ax.scatter(true_time-tshift + OFFSET - 2400000.5, true_time-est_time, c='k', s=0.5, label = "Phase Shift (BKJD Days)")
        if ii == 0:
            ax.set_ylabel('Epsilon (BKJD days)')
          
            pt.title(f"Epsilon Values for KIC {star_name}")
        ax.set_xlim(segments[ii][0]-tshift + OFFSET - 2400000.5 - margin, 
                    segments[ii][-1]-tshift + OFFSET - 2400000.5 + margin)
    
    positive_data = np.abs(true_time-est_time)
    m, b = np.polyfit(true_time-tshift + OFFSET - 2400000.5, positive_data, 1)
    x = true_time-tshift + OFFSET - 2400000.5
    pt.plot(x, x*m + b, color = 'r', label = "Linear Drift", linestyle = '--')
    fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
    """
    return true_time-est_time


def get_epsilon_value_tess_nav(star_name, sine_string_list):

    master_list_eps = []
    i = 0 
    for nameOfStar, sine_string in zip(star_name, sine_string_list):
        master_flux = []
        master_time = [] 
        #search_result = lk.search_lightcurve(f"{star_name}")
        #lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
        search_result = lk.search_tesscut(target=f"{nameOfStar}")
        search_result = search_result[0]
        print(search_result)
        tpf_collection = search_result.download_all(cutout_size=50)
        
        for l in tpf_collection:
            s = unpopular.Source(l.path, remove_bad=True)
            s.set_aperture(rowlims=[25, 26], collims=[25, 26])
            s.add_cpm_model(exclusion_size=5, n=64, predictor_method="similar_brightness")
            s.set_regs([0.1])
            s.holdout_fit_predict(k=100);

            aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux")
            aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction", weighting=None)
            #pt.plot(s.time, aperture_normalized_flux, ".", c="k", ms=8, label="Normalized Flux")
            #pt.plot(s.time, aperture_cpm_prediction, "-", lw=3, c="C3", alpha=0.8, label="CPM Prediction")
            #pt.xlabel("Time - 2457000 [Days]", fontsize=30)
            #pt.ylabel("Normalized Flux", fontsize=30)
            #pt.tick_params(labelsize=20)
            #pt.legend(fontsize=30)
            
            apt_detrended_flux = s.get_aperture_lc(data_type="cpm_subtracted_flux")
            min_val = np.percentile(apt_detrended_flux, 5)
            max_val = np.percentile(apt_detrended_flux,95)

            normalized_data = (2 * (apt_detrended_flux - min_val) / (max_val - min_val) - 1) 
            master_flux.extend(normalized_data)
            master_time.extend(s.time)
        
        lc = lk.TessLightCurve(time=s.time, flux=apt_detrended_flux) #lk.search_lightcurve(nameOfStar).download().remove_outliers(sigma = 5.0)
        #print(y.author)
    
    
    
        t = lc.time.value
        #sine_string = "0.0020 * np.sin(2* np.pi * (10.3376) * t + -0.2050) + 1 + 0.0017 * np.sin(np.pi * 2 * (12.4714) * t + -6.2832)"
        #sine_string = "0.0020 sin(2π(10.3376)t + -0.2050) + 1 + 0.0017 sin(2π(12.4714)t + -6.2832)"
        #0.0020  * np.sin(2 * np.pi * (10.3376) * t  + -0.2050) + 1 + 0.0017  * np.sin(2 * np.pi * (12.4714) * t  + -6.2832)
        sine_string = sine_string.replace('sin', 'np.sin')
        sine_string = sine_string.replace('2π', '2 * np.pi ')
        #sine_string = sine_string.replace('t', ' * t ')
        sine_string = sine_string.replace("f(t) = ", "")
        #model_profile = eval(sine_string)
        OFFSET = 2457000
        expected_cadence = 1426  # seconds

        
        def create_model_function(sine_string):
            """Create a callable function from the sine string"""
            def model(t, dt, *params):

                shifted_t = t + dt # + (OFFSET-2457000)
                print(sine_string)
                return eval(sine_string.replace('t', 'shifted_t'))
            return model

        
        #sine_string = "0.0020 * np.sin(2* np.pi * (10.3376) * t + -0.2050) + 1 + 0.0017 * np.sin(np.pi * 2 * (12.4714) * t + -6.2832)"
        profile_func = create_model_function(sine_string)

        #mask = (np.isfinite(lc.flux.value.unmasked))
        all_flux = lc.flux.value#[mask]
        all_time = lc.time.value#[mask]
        #print(all_time)
        true_time = []
        est_time = []
        t_step = 0.1  
        dt = 1.0  
        t_prev = -np.inf

        for t in all_time:
            if t - t_prev > t_step:
                t_prev = t
            else:
                continue

            
            start_bxjd = t_prev
            mask = (all_time >= start_bxjd) & (all_time <= (start_bxjd + dt))
            time = all_time[mask]
            flux = all_flux[mask]
            print(time)
            print(flux)
        
            if len(time) == 0:
                continue
            
            if abs(time[-1] - start_bxjd - dt) > expected_cadence/86400:
                continue
            if abs(time[0] - start_bxjd) > expected_cadence/86400:
                continue
            if np.any(np.diff(time) > expected_cadence/86400):
                continue
            
        
            t_zeroed = time - time[0]
            
            
            t0_list = np.arange(-4,4) * 0.1 + time[0] + np.random.normal(0.01, 0.05)
            t_est_list = []
            
            for t0 in t0_list:
                try:
                    popt, pcov = curve_fit(
                        profile_func,
                        xdata=t_zeroed,
                        ydata=flux,
                        p0=t0,
                        xtol=1e-12,
                        maxfev=1000
                    )
                    t_est_list.append(popt[0])
                except RuntimeError:
                    continue
            
            if len(t_est_list) > 0:
                t_est = t_est_list[np.argmin(np.abs(t_est_list-time[0]))]
                true_time.append(time[0])
                est_time.append(t_est)

    
        true_time = np.array(true_time)
        est_time = np.array(est_time)

        
        time_diff = np.diff(true_time)
        mask = time_diff > 3000
        gap_indices = np.where(mask)[0]
        segments = np.split(true_time, gap_indices+1)
        """
        # plotting
        tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5)/100)*100)
        margin = 0.5  # days
        
        fig_oc, axs = pt.subplots(1, len(segments), figsize=(8,3), sharey=True, 
                                gridspec_kw={'wspace': 0, 'hspace': 0},
                                width_ratios=[seg[-1]-seg[0] + margin*2 for seg in segments])

        if not isinstance(axs, np.ndarray):
            axs = [axs]

        for ii, ax in enumerate(axs):
            ax.scatter(true_time-tshift + OFFSET - 2400000.5, true_time-est_time, c='k', s=0.5, label = "Phase Shift (BKJD Days)")
            if ii == 0:
                ax.set_ylabel('Epsilon (BKJD days)')
            
                pt.title(f"Epsilon Values for KIC {star_name}")
            ax.set_xlim(segments[ii][0]-tshift + OFFSET - 2400000.5 - margin, 
                        segments[ii][-1]-tshift + OFFSET - 2400000.5 + margin)
        
        positive_data = np.abs(true_time-est_time)
        m, b = np.polyfit(true_time-tshift + OFFSET - 2400000.5, positive_data, 1)
        x = true_time-tshift + OFFSET - 2400000.5
        pt.plot(x, x*m + b, color = 'r', label = "Linear Drift", linestyle = '--')
        fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
        """
        eps = true_time-est_time
        half = int(len(eps)/2)
        eps_1_half = eps[:half]
        eps_2_half = eps[half:]
        print(np.average(eps))
        print(np.average(eps_1_half))
        print(np.average(eps_2_half))
        percentage = np.abs(np.average(eps_2_half)-np.average(eps_1_half))/((np.average(eps_1_half) + np.average(eps_2_half))/2)
        percentage *= 100
        master_list_eps.append({"Name": star_name[i], "average eps": np.average(eps), "1 Half": np.average(eps_1_half), "2 Half": np.average(eps_2_half), "Percent Diff": percentage})
        i += 1
    df = pd.DataFrame(master_list_eps)
    df.to_csv('KeplerStarsOutput_EPS_VALS_updated_tess_nav.csv', index=False)
    print("\nResults saved to KeplerStarsOutput")
        #master_list_eps.append({"KIC": star_name[i], "average eps": np.average(eps), "1 Half": np.average(eps_1_half), "2 Half": np.average(eps_2_half), "Percent Diff": percentage})
        #return true_time-est_time


def get_csv_epsilon_value(csv_file_path): 
    print("Running")
    try:
        df = pd.read_csv(csv_file_path)
        
        #if 'TIC_ID' not in df.columns:
        #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
        
        KIC_list = df['KIC'].dropna().astype(str).tolist()
        FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
        i = 0 
        master_list_eps = []
        while( i < len(KIC_list)):
            eps = get_epsilon_value(KIC_list[i], FUNCTION_list[i])
            if(i>0):
                pt.close('all')
            
            half = int(len(eps)/2)
            eps_1_half = eps[:half]
            eps_2_half = eps[half:]
            print(np.average(eps))
            print(np.average(eps_1_half))
            print(np.average(eps_2_half))
            master_list_eps.append({"KIC": KIC_list[i], "average eps": np.average(eps), "1 Half": np.average(eps_1_half), "2 Half": np.average(eps_2_half)})
            i += 1
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_EPS_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
    except Exception as e:
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_EPS_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
        print(f"Error loading TIC IDs: {e}")
        return []


def find_valid_segments(all_time, all_flux, dt=1.0, t_step=0.1, expected_cadence_days=1800/86400):
    segments = []
    i = 0
    N = len(all_time)
    while i < N:
        t0 = all_time[i]
        t1 = t0 + dt
        j = np.searchsorted(all_time, t1)

        time_chunk = all_time[i:j]
        flux_chunk = all_flux[i:j]

        if len(time_chunk) >= 10:
            diffs = np.diff(time_chunk)
            if np.all(diffs < 1.5 * expected_cadence_days):  # allow minor jitter
                segments.append((time_chunk, flux_chunk))

        i = np.searchsorted(all_time, t0 + t_step)

    return segments
                
def load_tic_ids_from_csv(csv_file_path):
    try:
        df = pd.read_csv(csv_file_path)
        
        #if 'TIC_ID' not in df.columns:
        #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
        
        tic_list = df['KIC'].dropna().astype(str).tolist()
        return tic_list
    
    except Exception as e:
        print(f"Error loading TIC IDs: {e}")
        return []


def seriesofstarsTest(listofstars):
    results = []
    try:
        for star in listofstars:
            print(f"KIC {star}")
            function, lc, composite_strings = getCompositeSine2_second_test(f"KIC {star}")
            if (function[0] == -10):
                continue
            flux = lc.flux.value
            print(flux)
            print(function)
            time = lc.time.value
            min_length = min(len(flux), len(function))
            flux = flux[:min_length]
            time = time[:min_length]
            function = function[:min_length]
            residuals = flux - function
            mse = np.sum((residuals)**2)/len(flux)
            print(f"MSE: {np.sum((residuals)**2)/len(flux)}")
            results.append({'KIC': star, 'MSE': mse, 'Composite Function': composite_strings})
    except Exception as e: 
        df = pd.DataFrame(results)
        df.to_csv('KeplerStarsOutput.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
        return results
    df = pd.DataFrame(results)
    df.to_csv('KeplerStarsOutput.csv', index=False)
    print("\nResults saved to KeplerStarsOutput")

def seriesofstarsTest_time_error(listofstars):
    results = []
    try:
        for star in listofstars:
            print(f"KIC {star}")
            time = identifyPeaksOfLightcurves(f"KIC {star}", 100)
            results.append({'KIC': star, 'Error': time})
    except RuntimeError as e: 
        df = pd.DataFrame(results)
        df.to_csv('KeplerStarsOutput_time_error.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
        return results
    df = pd.DataFrame(results)
    df.to_csv('KeplerStarsOutput_time_error.csv', index=False)
    print("\nResults saved to KeplerStarsOutput")

def CompareTelescopes(nameOfStar, sine_string): 
    lc = lk.search_lightcurve(f"TIC {nameOfStar}").download_all().stitch().remove_outliers(sigma = 1.0)
    signal = lc.flux.value
    t = lc.time.value
    sine_string = sine_string.replace('sin', 'np.sin')
    sine_string = sine_string.replace('2π', '2 * np.pi ')
    #sine_string = sine_string.replace('t', ' * t ')
    sine_string = sine_string.replace("f(t) = ", "")
    model = eval(sine_string)

    def normalize_to_minus_one_to_one(arr):
        min_val = np.min(arr)
        max_val = np.max(arr)
        normalized = 2 * (arr - min_val) / (max_val - min_val) -1
        return normalized
    model = normalize_to_minus_one_to_one(model)
    signal = normalize_to_minus_one_to_one(signal)
    #pt.plot(t, model)
    #pt.plot(t,signal)
    #pt.show()
    def spectral_goodness_of_fit(signal, model):
        """
        Computes the spectral residual and normalized R^2_FFT goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - signal: numpy array, the observed signal
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_FFT: float
        """
      
        n = 15 
        b = [1.0 / n] * n
        a = 1
        signal = lfilter(b, a, signal)
        
        S_f = np.abs(np.fft.rfft(signal))
        M_f = np.abs(np.fft.rfft(model))
        
       
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
       
        S_bar = np.mean(S_f)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_FFT = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_FFT
    
    spec_res, R2 = spectral_goodness_of_fit(signal, model)
    print(spec_res, R2)
    return spec_res, R2



def SpectralResiduals(nameOfStar, sine_string): 
    lc = lk.search_lightcurve(f"KIC {nameOfStar}").download_all().stitch().remove_outliers(sigma = 5.0)
    signal = lc.flux.value
    t = lc.time.value
    sine_string = sine_string.replace('sin', 'np.sin')
    sine_string = sine_string.replace('2π', '2 * np.pi ')
    #sine_string = sine_string.replace('t', ' * t ')
    sine_string = sine_string.replace("f(t) = ", "")
    model = eval(sine_string)
    """
    def normalize_to_minus_one_to_one(arr):
        min_val = np.min(arr)
        max_val = np.max(arr)
        normalized = 2 * (arr - min_val) / (max_val - min_val) -1
        return normalized
    """
    #model = normalize_to_minus_one_to_one(model)
    #signal = normalize_to_minus_one_to_one(signal)
    #pt.plot(t, model)
    #pt.plot(t,signal)
    #pt.show()
    def spectral_goodness_of_fit(signal, model):
        """
        Computes the spectral residual and normalized R^2_FFT goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - signal: numpy array, the observed signal
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_FFT: float
        """
      
        S_f = np.abs(np.fft.rfft(signal))
        M_f = np.abs(np.fft.rfft(model))
        
       
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
        
       
        S_bar = np.mean(S_f)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_FFT = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_FFT
    
    spec_res, R2 = spectral_goodness_of_fit(signal, model)
    print(spec_res, R2)
    return spec_res, R2
    
def SpectralResidualsCsvBased(csv_file_path): 
    print("Running")
    try:
        df = pd.read_csv(csv_file_path)
        
        #if 'TIC_ID' not in df.columns:
        #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
        
        KIC_list = df['KIC'].dropna().astype(str).tolist()
        FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
        i = 0 
        master_list_eps = []
        while( i < len(KIC_list)):
            spectral_resid, R2 = SpectralResiduals(KIC_list[i], FUNCTION_list[i])
            if(i>0):
                pt.close('all')
            
            #half = int(len(eps)/2)
            #eps_1_half = eps[:half]
            #eps_2_half = eps[half:]
            #print(np.average(eps))
            #print(np.average(eps_1_half))
            #print(np.average(eps_2_half))
            master_list_eps.append({"KIC": KIC_list[i], "Spectral Residuals": spectral_resid, "R2_FFT": R2})
            i += 1
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_Spectral_residual_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput_Spectral_residual_VALS")
    except Exception as e:
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_Spectral_residual_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput_Spectral_residual_VALS")
        print(f"Error loading TIC IDs: {e}")
        return []
    
def CompareTelescopesCsvBased(csv_file_path): 
    print("Running")
    try:
        df = pd.read_csv(csv_file_path)
        
        #if 'TIC_ID' not in df.columns:
        #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
        
        TIC_list = df['TIC'].dropna().astype(str).tolist()
        KIC_list = df['KIC'].dropna().astype(str).tolist()
        FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
        i = 0 
        master_list_eps = []
        while( i < len(KIC_list)):
            try:
                pt.close('all')
                search_result = lk.search_lightcurve(f"TIC {TIC_list[i]}", author='TASOC' )
                lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
                time = lc.time.value
                flux = lc.flux_corr.value
                #print(lc.time.unit)
                #lc.plot()
                pt.plot(time, flux)
                pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")
                pt.savefig(fr'C:\Users\ahmed\Downloads\data_images\picture_TIC_{TIC_list[i]}')
                #pt.xlabel(f'{lc.time.units}')
                #pt.ylabel(f'{lc.flux_corr.units}')
                
                #pt.show()
                #spec_res, R2 = CompareTelescopes(TIC_list[i], FUNCTION_list[i])
                #if(i>0):
                #    pt.close('all')
                
                #master_list_eps.append({"KIC": KIC_list[i], "TIC":TIC_list[i], "Spectral_Residuals": spec_res, "R2_FFT": R2})
                i += 1
            except: 
                i+= 1
                continue
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_TESS_Residual_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
    except Exception as e:
        df = pd.DataFrame(master_list_eps)
        df.to_csv('KeplerStarsOutput_TESS_Residual_VALS.csv', index=False)
        print("\nResults saved to KeplerStarsOutput")
        print(f"Error loading TIC IDs: {e}")
        return []
    
def plotMap():


    pt.style.use(['science', 'no-latex'])
    #pt.rcParams.update({'figure.dpi': '500'})

    with open(r"C:\Users\ahmed\Downloads\asu.tsv", 'r') as file:
        lines = file.readlines()

    kic_list = []
    ra_list = []
    de_list = []

    in_data_section = False

    for line in lines:
        line = line.strip()
        
       
        if line.startswith('#') or not line:
            in_data_section = False
            continue
        
        
        if line.startswith('_RAJ2000;_DEJ2000;KIC;RAJ2000;DEJ2000'):
            in_data_section = True
            continue
        
       
        if line.startswith('deg;deg; ;deg;deg'):
            continue
        
        
        if line.startswith('------------;------------;--------;----------;----------'):
            continue
        
        if in_data_section:
            parts = line.split(';')
            if len(parts) >= 5:
            
                ra = float(parts[3].strip())
                de = float(parts[4].strip())
                kic = int(parts[2].strip())
                
                ra_list.append(ra)
                de_list.append(de)
                kic_list.append(kic)


    kic_array = np.array(kic_list)
    ra_array = np.array(ra_list)
    de_array = np.array(de_list)
    df = pd.DataFrame({
    'KIC': kic_array,
    'Right Ascension': ra_array,
    'Declination': de_array})
    df.to_csv('star_data_locations.csv', index=False)

    
    print("KIC array:", kic_array[:5])
    print("RA array (degrees):", ra_array[:5])
    print("DE array (degrees):", de_array[:5])
    print(f"Total entries: {len(kic_array)}")

    #pt.figure(figsize=(12, 7))
    pt.scatter(ra_array, de_array, s=2, label='Analyzed Stars')
    pt.xlabel('Right Ascension (degrees)')
    pt.ylabel('Declination (degrees)')
    pt.title('Sky Position of Analyzed Stars')
    #plt.legend()
    pt.grid(True)
    pt.gca().invert_xaxis()
    pt.show()

    ra_shifted = np.where(ra_array > 180, ra_array - 360, ra_array)

    
    ra_rad = np.radians(ra_shifted)
    de_rad = np.radians(de_array)
    
    #pt.figure(figsize=(12, 7))
    pt.figure(figsize=(12, 7))
    ax = pt.subplot(111, projection="aitoff")
    pt.scatter(ra_rad, np.radians(de_array), s=1)
    pt.title('Analyzed Stars in Aitoff Projection')
    
    
    pt.grid(True)
    pt.xlabel('Right Ascension (degrees)')
    pt.ylabel('Declination (degrees)')
    pt.savefig("my_plot.png")
    pt.show()
    
def plotMap():
    pt.style.use(['science', 'no-latex'])

    with open(r"C:\Users\ahmed\Downloads\asu.tsv", 'r') as file:
        lines = file.readlines()

    kic_list = []
    ra_list = []
    de_list = []

    in_data_section = False

    for line in lines:
        line = line.strip()
        if line.startswith('#') or not line:
            in_data_section = False
            continue
        if line.startswith('_RAJ2000;_DEJ2000;KIC;RAJ2000;DEJ2000'):
            in_data_section = True
            continue
        if line.startswith('deg;deg; ;deg;deg') or line.startswith('------------'):
            continue

        if in_data_section:
            parts = line.split(';')
            if len(parts) >= 5:
                ra = float(parts[3].strip())
                de = float(parts[4].strip())
                kic = int(parts[2].strip())

                ra_list.append(ra)
                de_list.append(de)
                kic_list.append(kic)

    kic_array = np.array(kic_list)
    ra_array = np.array(ra_list)
    de_array = np.array(de_list)

    df = pd.DataFrame({
        'KIC': kic_array,
        'Right Ascension': ra_array,
        'Declination': de_array})
    #df.to_csv('star_data_locations.csv', index=False)

    print("KIC array:", kic_array[:5])
    print("RA array (degrees):", ra_array[:5])
    print("DE array (degrees):", de_array[:5])
    print(f"Total entries: {len(kic_array)}")

    pt.scatter(ra_array, de_array, s=2, label='Analyzed Stars')
    pt.xlabel('Right Ascension (degrees)')
    pt.ylabel('Declination (degrees)')
    pt.title('Sky Position of Analyzed Stars')
    pt.grid(True)
    pt.gca().invert_xaxis()
    pt.show()

    ra_shifted = np.where(ra_array > 180, ra_array - 360, ra_array)
    ra_rad = np.radians(ra_shifted)
    de_rad = np.radians(de_array)

    pt.figure(figsize=(12, 7))
    ax = pt.subplot(111, projection="aitoff")
    ax.scatter(ra_rad, de_rad, s=1)

    
    xticks_deg = np.arange(-180, 181, 60)
    yticks_deg = np.arange(-90, 91, 30)
    xtick_labels = [f"{int(t)}" if abs(t) != 180 else '' for t in xticks_deg]
    ytick_labels = [f"{int(t)}" if abs(t) != 90 else '' for t in yticks_deg]
    ax.set_xticks(np.radians(xticks_deg))
    ax.set_xticklabels(xtick_labels)
    ax.set_yticks(np.radians(yticks_deg))
    ax.set_yticklabels(ytick_labels)

    ax.set_title('Analyzed Stars in Aitoff Projection')
    ax.grid(True)
    pt.xlabel('Right Ascension (degrees)')
    pt.ylabel('Declination (degrees)')
    #pt.savefig("my_plot.png")
    pt.show()

def unpopular_clean_tess(csv_path):
    print("Running")
    #test_list = [271959957, 268160106,159302678, 239311449]
    #test_list = str(test_list)
    df = pd.read_csv(csv_path)
    
    #if 'TIC_ID' not in df.columns:
    #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
    
    TIC_list = df['TIC'].dropna().astype(str).tolist()[50:]
    KIC_list = df['KIC'].dropna().astype(str).tolist()[50:]
    FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()[50:]
    i = 0 
    master_lists_tess_pop = []
 
    def spectral_goodness_of_fit(signal, model):
        """
        Computes the spectral residual and normalized R^2_FFT goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - signal: numpy array, the observed signal
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_FFT: float
        """
        
        #n = 5  # the larger n is, the smoother curve will be
        #b = [1.0 / n] * n
        #a = 1
        #signal = lfilter(b, a, signal)
        #pt.plot(np.arange(len(signal)), signal)
        #pt.show()
        S_f = np.abs(np.fft.rfft(signal))
        pt.plot(np.arange(len(S_f)),S_f/np.max(S_f), color = 'red', label = 'signal')
        M_f = np.abs(np.fft.rfft(model))
        pt.plot(np.arange(len(M_f)),M_f/np.max(M_f), color = 'orange', label = 'model')
        def FilterPeaksfft(fft, scalar = 0.4):
            fft /= np.max(fft)
            peaks, _= find_peaks(fft, prominence=np.max(fft) * scalar)#height=[np.max(fft) * 0.55, np.max(fft) * 1.1])
            if(len(peaks) == 0):
                peaks = [np.argmax(fft)]
            #pt.figure(figsize=(10, 6))
            #pt.plot(pg.frequency, pg.power, label='Periodogram')
            peak_amps = fft[peaks]
            filtered_fft = np.zeros(len(fft))
            #peak_index = peaks
            for peak_index in peaks:
                if peak_index < 100:
                    continue
                filtered_fft[peak_index - 25: peak_index+25] = fft[peak_index-25:peak_index+25]
            #pt.plot(np.arange(len(filtered_fft)), filtered_fft)
            #filtered_fft = np.fft.ifft(filtered_fft)
            #print(peak_amps)
            #pt.plot(np.arange(len(fft)), fft, color = 'green')
            #pt.scatter(peaks, peak_amps, color = 'red')
            
            #pt.show()
            
            #print(pg.frequency[filtered_peaks])
            #print(pg.power[filtered_peaks])
            return filtered_fft
        S_f = FilterPeaksfft(S_f)
        M_f = FilterPeaksfft(M_f, scalar = 0.2)
        signal = np.fft.ifft(S_f)

        pt.plot(np.arange(len(M_f)), S_f, label = 'Signal_clean')
        pt.plot(np.arange(len(M_f)), M_f, label = 'Model_clean')
        pt.title(f"KIC {KIC_list[i]} | TIC {TIC_list[i]}")
        pt.legend()
        #pt.show()
        pt.savefig(fr"C:\Users\ahmed\Downloads\higher try\pic_TIC_{TIC_list[i]}.png")
        pt.close()
        
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
       
        S_bar = np.mean(signal)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_FFT = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_FFT
    
    while( i < len(KIC_list)):
        try:
            master_flux = []
            master_time = []
            sine_string = FUNCTION_list[i]
            pt.close('all')
            #if TIC_list[i] not in test_list:
            #    i+= 1
            #    continue
            search_result = lk.search_tesscut(target=f"TIC{TIC_list[i]}", sector  = [14,15] )
            
            print(search_result)
            tpf_collection = search_result.download_all(cutout_size=50)
            
            for l in tpf_collection:
                s = unpopular.Source(l.path, remove_bad=True)
                s.set_aperture(rowlims=[25, 26], collims=[25, 26])
                s.add_cpm_model(exclusion_size=5, n=64, predictor_method="similar_brightness")
                s.set_regs([0.1])
                s.holdout_fit_predict(k=100);

                aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux")
                aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction", weighting=None)
                #pt.plot(s.time, aperture_normalized_flux, ".", c="k", ms=8, label="Normalized Flux")
                #pt.plot(s.time, aperture_cpm_prediction, "-", lw=3, c="C3", alpha=0.8, label="CPM Prediction")
                #pt.xlabel("Time - 2457000 [Days]", fontsize=30)
                #pt.ylabel("Normalized Flux", fontsize=30)
                #pt.tick_params(labelsize=20)
                #pt.legend(fontsize=30)
                
                apt_detrended_flux = s.get_aperture_lc(data_type="cpm_subtracted_flux")
                min_val = np.percentile(apt_detrended_flux, 5)
                max_val = np.percentile(apt_detrended_flux,95)

                normalized_data = (2 * (apt_detrended_flux - min_val) / (max_val - min_val) - 1) 
                master_flux.extend(normalized_data)
                master_time.extend(s.time)
            t = np.array(master_time)
            sine_string = sine_string.replace('sin', 'np.sin')
            #sine_string = sine_string.replace('2π', '2 * np.pi ')
            sine_string = sine_string.replace('2??', '2 * np.pi' )
            #sine_string = sine_string.replace('t', ' * t ')
            sine_string = sine_string.replace("f(t) = ", "")
            print(sine_string)
            model = eval(sine_string)
            model = 2 * (model - np.min(model)) / (np.max(model) - np.min(model)) - 1
            #model *= np.max(master_flux)
            spec, R2  = spectral_goodness_of_fit(master_flux, model )
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2FFT": R2})
            #pt.plot(master_time, master_flux, "k-")
            pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")
            #pt.savefig(fr"C:\Users\ahmed\Downloads\Filtering Results\unpop_pic_{TIC_list[i]}.png")
            i+=1

        except Exception as the_exception:
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": 0,"R2FFT": -1})
            print('did not find')
            print(the_exception)
            df = pd.DataFrame(master_lists_tess_pop)
            df.to_csv(r"C:\Users\ahmed\Downloads\higher try\results.csv")
            i+=1
            continue


        

    df = pd.DataFrame(master_lists_tess_pop)
    df.to_csv(fr"C:\Users\ahmed\Downloads\higher try\results.csv")


def tess_clean_MAST(csv_path):
    print("Running")

    df = pd.read_csv(csv_path)
    
    #if 'TIC_ID' not in df.columns:
    #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
    
    TIC_list = df['TIC'].dropna().astype(str).tolist()[:10]
    KIC_list = df['KIC'].dropna().astype(str).tolist()[:10]
    FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()[:10]
    i = 0 
    master_lists_tess_pop = []
 
    def spectral_goodness_of_fit(signal, model):
        """
        Computes the spectral residual and normalized R^2_FFT goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - signal: numpy array, the observed signal
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_FFT: float
        """
       
        #n = 15  # the larger n is, the smoother curve will be
        #b = [1.0 / n] * n
        #a = 1
        #signal = lfilter(b, a, signal)
        
        S_f = np.abs(np.fft.rfft(signal))
        M_f = np.abs(np.fft.rfft(model))
        
        
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
        
       
        S_bar = np.mean(S_f)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_FFT = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_FFT
    
    while( i < len(KIC_list)):
        try:
            #master_flux = []
            #master_time = []
            sine_string = FUNCTION_list[i]
            pt.close('all')
            

            search_result = lk.search_lightcurve(f"TIC {TIC_list[i]}", author='TASOC' )
            lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
            t = lc.time.value
            signal = lc.flux_corr.value

            signal = 2 * (signal - np.min(signal)) / (np.max(signal) - np.min(signal)) - 1

            t = np.array(t)
            sine_string = sine_string.replace('sin', 'np.sin')
            sine_string = sine_string.replace('2π', '2 * np.pi ')
            #sine_string = sine_string.replace('t', ' * t ')
            sine_string = sine_string.replace("f(t) = ", "")
            print(sine_string)
            model = eval(sine_string)
            model = 2 * (model - np.min(model)) / (np.max(model) - np.min(model)) - 1
            spec, R2  = spectral_goodness_of_fit(signal, model)
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2FFT": R2})
            pt.plot(t, signal, "k-")
            pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")
            pt.savefig(fr"C:\Users\ahmed\research_delta\ResearchPython\folder_tess_MAST\NATIVE_pic_{TIC_list[i]}.png")
            i+=1

        except:
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": 0,"R2FFT": -100000})
            print('did not find')
            i+=1
            continue


        

    df = pd.DataFrame(master_lists_tess_pop)
    df.to_csv(fr"C:\Users\ahmed\research_delta\ResearchPython\folder_tess_MAST\results.csv")

# only for plotting  
def unpopular_clean_tess_plotting(csv_path):
    print("Running")
    df = pd.read_csv(csv_path)

    TIC_list = df['TIC'].dropna().astype(str).tolist()
    KIC_list = df['KIC'].dropna().astype(str).tolist()
    FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
    i = 0 
    master_lists_tess_pop = []
    
    def spectral_goodness_of_fit(signal, model, master_time):
        dt = master_time[1] - master_time[0]
        S_f = np.abs(np.fft.rfft(signal))
        x_axis = np.fft.rfftfreq(n=2*len(S_f)-1, d = dt)
        pt.plot(x_axis,S_f/np.max(S_f), color = 'red', label = 'TESS Light Curve')
        M_f = np.abs(np.fft.rfft(model))
        pt.plot(x_axis,M_f/np.max(M_f), color = 'orange', label = 'Model')
        def FilterPeaksfft(fft, scalar = 0.4):
            fft /= np.max(fft)
            peaks, _= find_peaks(fft, prominence=np.max(fft) * scalar)#height=[np.max(fft) * 0.55, np.max(fft) * 1.1])
            if(len(peaks) == 0):
                peaks = [np.argmax(fft)]
            #pt.figure(figsize=(10, 6))
            #pt.plot(pg.frequency, pg.power, label='Periodogram')
            peak_amps = fft[peaks]
            filtered_fft = np.zeros(len(fft))
            #peak_index = peaks
            for peak_index in peaks:
                if peak_index < 100:
                    continue
                filtered_fft[peak_index - 25: peak_index+25] = fft[peak_index-25:peak_index+25]
            #pt.plot(np.arange(len(filtered_fft)), filtered_fft)
            #filtered_fft = np.fft.ifft(filtered_fft)
            #print(peak_amps)
            #pt.plot(np.arange(len(fft)), fft, color = 'green')
            #pt.scatter(peaks, peak_amps, color = 'red')
            
            #pt.show()
            
            #print(pg.frequency[filtered_peaks])
            #print(pg.power[filtered_peaks])
            return filtered_fft
        S_f = FilterPeaksfft(S_f)
        M_f = FilterPeaksfft(M_f, scalar = 0.2)
        signal = np.fft.ifft(S_f)

        pt.plot(x_axis, S_f, label = 'Filtered TESS Light Curve')
        pt.plot(x_axis, M_f, label = 'Filtered Model')
        pt.title(f"FFT for KIC {KIC_list[i]} and TIC {TIC_list[i]}")
        pt.xlabel('Frequency (cycles/BKJD)')
        pt.ylabel('Normalized Power')
        pt.legend()
        #pt.show()
        #pt.savefig(fr"C:\Users\ahmed\Downloads\higher try\pic_TIC_{TIC_list[i]}.png")
        pt.show()
        
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
       
        
        S_bar = np.mean(signal)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_FFT = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_FFT
    
    while( i < len(KIC_list)):
        try:
            master_flux = []
            master_time = []
            sine_string = FUNCTION_list[i]
            pt.close('all')
            #if TIC_list[i] not in test_list:
            #    i+= 1
            #    continue
            search_result = lk.search_tesscut(target=f"TIC {TIC_list[i]}", sector  = [14,15] )
            
            print(search_result)
            tpf_collection = search_result.download_all(cutout_size=50)
            
            for l in tpf_collection:
                s = unpopular.Source(l.path, remove_bad=True)
                s.set_aperture(rowlims=[25, 26], collims=[25, 26])
                s.add_cpm_model(exclusion_size=5, n=64, predictor_method="similar_brightness")
                s.set_regs([0.1])
                s.holdout_fit_predict(k=100);

                aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux")
                aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction", weighting=None)
                #pt.plot(s.time, aperture_normalized_flux, ".", c="k", ms=8, label="Normalized Flux")
                #pt.plot(s.time, aperture_cpm_prediction, "-", lw=3, c="C3", alpha=0.8, label="CPM Prediction")
                #pt.xlabel("Time - 2457000 [Days]", fontsize=30)
                #pt.ylabel("Normalized Flux", fontsize=30)
                #pt.tick_params(labelsize=20)
                #pt.legend(fontsize=30)
                
                apt_detrended_flux = s.get_aperture_lc(data_type="cpm_subtracted_flux")
                min_val = np.percentile(apt_detrended_flux, 5)
                max_val = np.percentile(apt_detrended_flux,95)

                normalized_data = (2 * (apt_detrended_flux - min_val) / (max_val - min_val) - 1) 
                master_flux.extend(normalized_data)
                master_time.extend(s.time)
            t = np.array(master_time)
            sine_string = sine_string.replace('sin', 'np.sin')
            #sine_string = sine_string.replace('2π', '2 * np.pi ')
            sine_string = sine_string.replace('2??', '2 * np.pi' )
            #sine_string = sine_string.replace('t', ' * t ')
            sine_string = sine_string.replace("f(t) = ", "")
            print(sine_string)
            model = eval(sine_string)
            model = 2 * (model - np.min(model)) / (np.max(model) - np.min(model)) - 1
            #model *= np.max(master_flux)
            spec, R2  = spectral_goodness_of_fit(master_flux, model, master_time )
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2FFT": R2})
            #pt.plot(master_time, master_flux, "k-")
            pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")
            #pt.savefig(fr"C:\Users\ahmed\Downloads\Filtering Results\unpop_pic_{TIC_list[i]}.png")
            i+=1

        except Exception as the_exception:
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": 0,"R2FFT": -1})
            print('did not find')
            print(the_exception)
            #df = pd.DataFrame(master_lists_tess_pop)
            #df.to_csv(r"C:\Users\ahmed\Downloads\higher try\results.csv")
            i+=1
            continue


        

    #df = pd.DataFrame(master_lists_tess_pop)
    #df.to_csv(fr"C:\Users\ahmed\Downloads\higher try\results.csv")




#print(find_max_frequency("f(t) = 0.0047 * sin(2π * 1.0983 * t + 2.7960) + 0.7343 + 0.0005 * sin(2π * 1.5464 * t + 3.9936) + 0.0798 + 0.0012 * sin(2π * 2.0353 * t + 3.7956) + 0.1853"))
#print(percenterror(-216, -258 ))

#getPeriodogramData('')
#print(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\tic_ids.csv"))

#fig, a = pt.subplots(1)
#getLightCurveData("KIC 3733346")
#identifyPeaks('X Caeli')
#identifyPeaks('X Caeli')

#pt.rc('legend', labelsize= 12)

#pt.style.use(['science', 'no-latex'])
#pt.rcParams.update({'figure.dpi': '300'})
a = ['KIC 12602250', 'KIC 9700322','KIC 4048494', 'KIC 6951642', 'KIC 8623953', 'KIC 8197761']
b = ['KIC 3429637', 'KIC 10451090', 'KIC 2987660']
c = ['KIC 12602250' , 'KIC 9700322' , 'KIC 8197761' , 'KIC 8623953' ,  'KIC 6382916' ,'KIC 3429637']
d = ['2987660' , '10451090' , '8197761' , '8623953' ,'3429637']
e = ['2987660' ,'3429637']

#KIC 8197761!!!!!!!!' 'V593 Lyr',
#getChiSquaredReduced('BO Lyn')
#print(getChiSquared('KIC 8197761'))
#('KIC 3429637', 0) 
#12602250
#g, h = compGetPeriodogramData('KIC 2168333') ###
#g.plot()
somestars = ['9652324', '9653684', '9835555', '9653838', '9835690', '9593010', '9531207', '9654046', '9654088', '9531257', '9469972', '9531319', '9654221', '9593346', '9836103', '9593399', '9775887', '9531736', '9531803', '9593837', '9593892', '9896552', '9593927', '9532154', '9715991', '9532219', '9896704', '9896727', '9594189', '9655316', '9837085', '9471324', '9716440', '9897434', '9716778', '9594890', '9656027', '9897710', '9777444', '9717148', '9838223', '9717468', '9717684', '9717781', '9595900', '9717970', '9596313', '9778648', '9839247', '9473636', '9473646', '9779140', '9718903', '9899540', '9473963', '9839899', '9779523', '9719477']
somestars = ["9653684", "9469972", "9531319", "9775887", "9593837", "9896552", "9715991", "9532219", "9594189","9655316", "9716440", "9594890", "9897710", "9777444","9717148", "9717468", "9595900", "9717970", "9596313", "9779523"]

#"9717781" TESST SEPERATE
#seriesofstarsTest(somestars)
#results = []
#print(len(somestars))
#for i in load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\tic_ids.csv"): 
 #   print(i)
  #  g, h = compGetPeriodogramData(f'KIC {i}') 
   # g.plot()
    #pt.show()
#print(results)
#get_epsilon_value()
#eps = get_epsilon_value("3429637", "f(t) = 0.0022 * sin(2π * 10.3376 * t + -0.2375) + 0.4865 + 0.0005 * sin(2π * 10.9363 * t + -6.2832) + 0.1066 + 0.0018 * sin(2π * 12.4714 * t + -6.2832) + 0.4069")

#eps = get_epsilon_value("12268220", "f(t) = 0.0006 * sin(2π * 1.3580 * t + -0.4060) + 0.2128 + 0.0005 * sin(2π * 1.8051 * t + 6.2832) + 0.1782 + 0.0004 * sin(2π * 2.2606 * t + 1.5659) + 0.1511 + 0.0003 * sin(2π * 2.7156 * t + -1.9018) + 0.1190 + 0.0003 * sin(2π * 3.1640 * t + -1.6757) + 0.0954 + 0.0002 * sin(2π * 22.1534 * t + -0.8883) + 0.0912 + 0.0004 * sin(2π * 23.6299 * t + 1.6420) + 0.1505")
#stars_named = ['GN And', 'V0340 And', 'KW Aur', 'iot Boo', 'kap 2 Boo', 'AO CVn', 'bet Cas', 'V0701 CrA', 'del Del', 'LM Hya']
#stars_named = ['del Del', 'LM Hya']
#for val in stars_named:
#   #try: 
#    plotsidebysideactual(val)
#   #except Exception as e: 
#    #    print(e)
    
#print(getCompositeSine2_second_test("12602250"))
#print(getCompositeSine2_second_test("12268220"))
#half = int(len(eps)/2)
#eps_1_half = eps[:half]
#eps_2_half = eps[half:]
#print(np.average(eps))
#print(np.average(eps_1_half))
#print(np.average(eps_2_half))
#plotsidebysideactual("KIC 8197761")
#unpopular_clean_tess_plotting(r"C:\Users\ahmed\research_delta\KeplerStarsOutput_TIC_enabled_with_TIC_new.csv")
#seriesofstarsTest_time_error(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\KeplerStarsOutput_2_timeerror.csv"))
#identifyPeaksOfLightcurves_manual('KIC 3123138', 0)
#guessLegacy('KIC 4048494',0) 
#print(getMeanSquaredResidual('KIC 7548479'))
#guessIterative('KIC 3429637',0)
#print(guess("X Caeli", 1))
#3429637
#CompareTelescopes("63122689", "f(t) = 0.1970 * sin(2π * 3.6556 * t + 1.9644) + 0.8538 + 0.0360 * sin(2π * 7.3117 * t + -6.2832) + 0.1561")
#get_csv_epsilon_value(r"KeplerStarsOutput.csv")
#identifyPeaksPowerComp('3429637')
#print(pt.rcParams.keys())
  # Set y-tick label font size
#plotMap()
#print(SpectralResiduals("12268220", "f(t) = 0.0006 * sin(2π * 1.3580 * t + -0.4060) + 0.2128 + 0.0005 * sin(2π * 1.8051 * t + 6.2832) + 0.1782 + 0.0004 * sin(2π * 2.2606 * t + 1.5659) + 0.1511 + 0.0003 * sin(2π * 2.7156 * t + -1.9018) + 0.1190 + 0.0003 * sin(2π * 3.1640 * t + -1.6757) + 0.0954 + 0.0002 * sin(2π * 22.1534 * t + -0.8883) + 0.0912 + 0.0004 * sin(2π * 23.6299 * t + 1.6420) + 0.1505"))

#lc = lk.search_lightcurve("TIC 137817459").download_all().stitch().remove_outliers(sigma = 5.0)
#lc_kep = lk.search_lightcurve("KIC 2581626").download_all().stitch().remove_outliers(sigma =5.0)
#pt.figure(figsize=(10, 6))

#lc = lc_kep
#t = lc.time.value
#lc.plot() 
#lc_kep.plot()
#pul_1 = 0.07190 * np.sin(2*np.pi*(18.139700)*t + -2.440400) 
#pul_2 = 0.015700*np.sin(2*np.pi*(23.498500)*t + -0.000015)
#pul_3 = 0.007800 * np.sin(2*np.pi*(5.358700)*t - 1.669385)  + 1.004000
#pul_4 = pul_1 + pul_2 + pul_3 
#lc.plot()
#print(lc.time.value[0])
#pt.plot(lc.time.value, pul_4, label = "Pulsation Mode 1")
#t.plot(lc.time.value, lc.flux.value, label = "Light Curve")
#pt.plot(lc.time.value, pul_2, label = "Pulsation Mode 2")
#pt.plot(lc.time.value, pul_3, label = "Pulsation Mode 3")
#pt.title("Pulsation Modes for KIC 3429637")
#pt.xlabel("Time -2454833 [BKJD Days]")
#pt.ylabel("Normalized Flux")
#pt.legend()
#getPeriodogramData('KIC 3429637')
#identifyPeaks('KIC 3429637')
#seriesofstarsTest(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\tic_ids.csv"))
pt.show()


"""
81976
f(t) = 0.0047 * sin(2π * 1.0983 * t + 2.7960) + 0.7343 + 0.0005 * sin(2π * 1.5464 * t + 3.9936) + 0.0798 + 0.0012 * sin(2π * 2.0353 * t + 3.7956) + 0.1853')
-0.01489508114251666 avg
-0.01961746417338575
-0.010174312049798923
0.9999577448689809

12602250
f(t) = 0.0004 * sin(2π * 11.6214 * t + 1.8740) + 0.6758 + 0.0002 * sin(2π * 14.9794 * t + 3.3672) + 0.3244')
0.009388682743293237
0.0017640377816918397
0.00011398020267873233
0.9998491956880848

12268220
f(t) = 0.0006 * sin(2π * 1.3580 * t + -0.4060) + 0.2128 + 0.0005 * sin(2π * 1.8051 * t + 6.2832) + 0.1782 + 0.0004 * sin(2π * 2.2606 * t + 1.5659) + 0.1511 + 0.0003 * sin(2π * 2.7156 * t + -1.9018) + 0.1190 + 0.0003 * sin(2π * 3.1640 * t + -1.6757) + 0.0954 + 0.0002 * sin(2π * 22.1534 * t + -0.8883) + 0.0912 + 0.0004 * sin(2π * 23.6299 * t + 1.6420) + 0.1505')
0.0030193749036970486
0.0050080278373302225
0.0010307219700638748
0.9997062116709856

"""
"""

11027806 OMIT



KIC 2297728
Could not get data, lightcurve has corrupted files???? No idea. Only experienced this once prior. 

KIC 9353572
MSE: 1.3938312779313455e-06
RMSE: 0.515125254394959
f(t) = 0.0003 * sin(2π * 10.8822 * t + 0.2737) + 0.1476 + 0.0010 * sin(2π * 13.3928 * t + 1.1594) + 0.5553 + 0.0001 * sin(2π * 13.8161 * t + -0.0225) + 0.0651 + 0.0004 * sin(2π * 15.7575 * t + -0.8967) + 0.2321
Epsilon: 0.006914152545426193


KIC 2304168
MSE: 0.00013053541098047208
RMSE: 0.7005639971786805
f(t) = 0.0080 * sin(2π * 8.1184 * t + 1.6094) + 0.5430 + 0.0032 * sin(2π * 8.5925 * t + 0.1283) + 0.2148 + 0.0026 * sin(2π * 10.4877 * t + 2.6149) + 0.1761 + 0.0010 * sin(2π * 11.0852 * t + -0.3035) + 0.0665
Epsilon: 0.06964311105059576

KIC 3123138
MSE: 1.5358083863541634e-06
RMSE: 0.7671053729865529
f(t) = 0.0002 * sin(2π * 1.0098 * t + -0.9006) + 0.1042 + 0.0001 * sin(2π * 2.0404 * t + -2.5190) + 0.1031 + 0.0001 * sin(2π * 3.0706 * t + -3.3178) + 0.0654 + 0.0002 * sin(2π * 9.7877 * t + 1.3097) + 0.1154 + 0.0009 * sin(2π * 15.1035 * t + -1.3249) + 0.6119
Epsilon: 0.010544980738373127 --> Inspecting chart, it seems like it fluctuates between 0 error to a constant error depending on the region within the light curve

"""



"""
stars_processed =  ['GN And', 'V0340 And', 'KW Aur', 'iot Boo', 'kap 2 Boo', 'AO CVn', 'bet Cas', 'del Del', 'LM Hya']
functions_processed = ["f(t) = 0.0063 * sin(2π * 14.4282 * t + 0.7797) + -0.0002", "f(t) = 0.0051 * sin(2π * 22.2584 * t + 6.2832) + 0.0004", "f(t) = 0.0010 * sin(2π * 11.3495 * t + -6.2832) + -0.0000", "f(t) = 0.0002 * sin(2π * 1.0987 * t + 6.2832) + 0.0002 + 0.0001 * sin(2π * 11.2529 * t + -6.2832) + 0.0002 + 0.0001 * sin(2π * 18.3035 * t + -6.2832) + 0.0002", "f(t) = 0.0002 * sin(2π * 1.1884 * t + 6.2832) + 0.0001 + 0.0001 * sin(2π * 1.5484 * t + -6.2832) + 0.0001 + 0.0001 * sin(2π * 1.9670 * t + -6.2832) + 0.0001 + 0.0001 * sin(2π * 7.0882 * t + -6.2832) + 0.0001 + 0.0001 * sin(2π * 14.6102 * t + -6.2832) + 0.0001 + 0.0002 * sin(2π * 15.4326 * t + -6.2832) + 0.0001", "f(t) = 0.0005 * sin(2π * 1.5106 * t + -6.2832) + -0.0000 + 0.0004 * sin(2π * 8.2471 * t + -6.2832) + -0.0000", "f(t) = 0.0011 * sin(2π * 9.8970 * t + -6.2832) + 0.0002", "f(t) = 0.0001 * sin(2π * 1.0756 * t + -6.2832) + -0.0000 + 0.0001 * sin(2π * 1.5539 * t + 6.2832) + -0.0000 + 0.0001 * sin(2π * 2.3647 * t + -6.2832) + -0.0000 + 0.0001 * sin(2π * 6.5413 * t + -6.2832) + -0.0000 + 0.0001 * sin(2π * 12.7654 * t + -6.2832) + -0.0000", "f(t) = 0.0002 * sin(2π * 1.2283 * t + 6.2832) + 0.0000 + 0.0002 * sin(2π * 21.7072 * t + -6.2832) + 0.0000"]

get_epsilon_value_tess_nav(stars_processed, functions_processed)
#print(len(stars_processed))
#print(len(functions_processed))
"""