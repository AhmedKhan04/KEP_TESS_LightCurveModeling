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
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import interp1d

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
    x = lk.search_targetpixelfile(nameOfStar).download().to_lightcurve()
    y = lk.search_lightcurve(nameOfStar).download_all().stitch().remove_outliers(sigma = 5.0)
    ## quarter term ^^^  y = lk.search_lightcurve(nameOfStar, quarter=(6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    z = x.to_periodogram()
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

def identifyPeaks(nameOfStar, lowerscalar = 0.1):
    pg, lightc = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * lowerscalar, max_power * 1.1])
    pt.figure(figsize=(10, 6))
    #pt.plot(pg.frequency, pg.power, label='Periodogram')
    x = pg.frequency[peaks]
    #for i in x:
     #   print(i.value)
    y = pg.power[peaks]
    filtered_peaks = []
    # Iterate through the peaks and check for closeness
    for i in range(len(x)):
        if len(filtered_peaks) == 0:
            if(x[i].value >= 1):
                filtered_peaks.append(peaks[i])
        else:
            # current peak is close to the last added peak?
            if np.abs(x[i].value - pg.frequency[filtered_peaks[-1]].value) <= 0.3:
                if(y[i].value > pg.power[filtered_peaks[-1]].value):
                    filtered_peaks[-1]= peaks[i]
            else: 
                filtered_peaks.append(peaks[i])
    if(len(pg.frequency[filtered_peaks]) > 10):
        return -1, 0        
    #pt.scatter(pg.frequency[filtered_peaks], pg.power[filtered_peaks], color='red', zorder=5, label='Local Maxima')
    #pt.xlabel('Frequency (cycles/day)')
    #pt.ylabel('Power')
    #pt.title('Periodogram with Local Maxima: '+ nameOfStar)
    #pt.legend()
    #pt.show()
    return(pg.frequency[filtered_peaks], lightc, pg.power[filtered_peaks])

def identifyPeaksPowerComp(nameOfStar):
    pg, ltcurves = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * 0.1, max_power * 1.1])
    #pt.figure(figsize=(10, 6))
    #pt.plot(pg.frequency, pg.power, label='Periodogram')
    x = pg.frequency[peaks]
    #for i in x:
     #   print(i.value)
    y = pg.power[peaks]
    filtered_peaks = []
    # Iterate through the peaks and check for closeness
    for i in range(len(x)):
        if len(filtered_peaks) == 0:
            if(x[i].value >= 1):
                filtered_peaks.append(peaks[i])
        else:
            # current peak is close to the last added peak?
            if np.abs(x[i].value - pg.frequency[filtered_peaks[-1]].value) <= 0.3:
                if(y[i].value > pg.power[filtered_peaks[-1]].value):
                    filtered_peaks[-1]= peaks[i]
            else: 
                filtered_peaks.append(peaks[i])
    if(len(pg.frequency[filtered_peaks]) > 10):
        return -1, 0
    #pt.scatter(pg.frequency[filtered_peaks], pg.power[filtered_peaks], color='red', zorder=5, label='Local Maxima')
    #pt.xlabel('Frequency (cycles/day)')
    #pt.ylabel('Power')
    #pt.title('Periodogram with Local Maxima: '+ nameOfStar)
    #pt.legend()
    #pt.show()
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
        #WORK AND FIX THIS PART

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
        #Neutral 

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
            
            # Prepare initial guesses and bounds
            ig_low = [low_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            ig_high = [high_amplitude_guess_scale * amplitude_guess, phase_guess, frequency_guess, offset_guess]
            
            bounds_low = ([low_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], 
                        [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.percentile(flux,95)])
            bounds_high = ([high_amplitude_scale*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.percentile(flux,5)], 
                        [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.percentile(flux,95)])
            
            # Try all three fits (low, current, high)
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
            
            # Compare all three options
            options = [
                ('low', fit_low_MSE, params_low, low_amplitude_scale),
                ('current', current_MSE, params, amplitude_scale),
                ('high', fit_high_MSE, params_high, high_amplitude_scale)
            ]
            
            # Find the best option
            best_option = min(options, key=lambda x: x[1])
            
            if best_option[0] == 'current':
                bestFitAchieved = True
                print("-1")
            else:
                # Update to the better option
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
        #Neutral 

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
        if(powerOfPeaks == -1):
            return [-10], 0
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

    search_result = lk.search_lightcurve(star_name, quarter=(6,7,8))
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
    OFFSET = 2454833
    expected_cadence = 1800  # seconds

    # Your model definition
    def create_model_function(sine_string):
        """Create a callable function from the sine string"""
        def model(t, dt, *params):
            # dt is the time offset parameter we're fitting
            # params could be used if you want to make amplitudes/frequencies variable
            shifted_t = t + dt + (OFFSET-2457000)
            print(sine_string)
            return eval(sine_string.replace('t', 'shifted_t'))
        return model

    # Define your model
    #sine_string = "0.0020 * np.sin(2* np.pi * (10.3376) * t + -0.2050) + 1 + 0.0017 * np.sin(np.pi * 2 * (12.4714) * t + -6.2832)"
    profile_func = create_model_function(sine_string)

    mask = (np.isfinite(lc.flux.value.unmasked))
    all_flux = lc.flux.value[mask]
    all_time = lc.time.value[mask]

    true_time = []
    est_time = []
    t_step = 0.1  # resolution of the Allan variance plot
    dt = 1.0  # days, duration of measurement for curve fit
    t_prev = -np.inf

    for t in all_time:
        if t - t_prev > t_step:
            t_prev = t
        else:
            continue

        # grab a segment of actual light curve data
        start_bxjd = t_prev
        mask = (all_time >= start_bxjd) & (all_time <= (start_bxjd + dt))
        time = all_time[mask]
        flux = all_flux[mask]

        # check data quality
        if len(time) == 0:
            continue
        if abs(time[-1] - start_bxjd - dt) > expected_cadence/86400:
            continue
        if abs(time[0] - start_bxjd) > expected_cadence/86400:
            continue
        if np.any(np.diff(time) > expected_cadence/86400):
            continue

        t_zeroed = time - time[0]
        
        # Initial guesses near the start time
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

    # Rest of your plotting code remains the same...
    true_time = np.array(true_time)
    est_time = np.array(est_time)

    # detect gaps and create segments
    time_diff = np.diff(true_time)
    mask = time_diff > 3000
    gap_indices = np.where(mask)[0]
    segments = np.split(true_time, gap_indices+1)

    # plotting
    tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5)/100)*100)
    margin = 0.5  # days

    fig_oc, axs = pt.subplots(1, len(segments), figsize=(8,3), sharey=True, 
                            gridspec_kw={'wspace': 0, 'hspace': 0},
                            width_ratios=[seg[-1]-seg[0] + margin*2 for seg in segments])

    if not isinstance(axs, np.ndarray):
        axs = [axs]

    for ii, ax in enumerate(axs):
        ax.scatter(true_time-tshift + OFFSET - 2400000.5, true_time-est_time, c='k', s=1)
        if ii == 0:
            ax.set_ylabel('O-C (days)', fontsize=11)
        ax.set_xlim(segments[ii][0]-tshift + OFFSET - 2400000.5 - margin, 
                    segments[ii][-1]-tshift + OFFSET - 2400000.5 + margin)

    fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
    pt.show()

    """
    expected_cadence = 1800 # seconds
    dt = 10  # days, duration of each segment for fitting
    t_step = 1  # days, step size for segments

    # Prepare cleaned time and flux arrays
    mask = np.isfinite(lc.flux.value)
    all_time = lc.time.value[mask]
    all_flux = lc.flux.value[mask]
    print(len(all_time))
    # Your precomputed model arrays (example)
    model_time = np.linspace(0, dt, len(model_profile))
    model_flux = model_profile

    model_interp = interp1d(model_time, model_flux, bounds_error=False, fill_value='extrapolate')

    true_time = []
    est_time = []

    # Create regularly spaced segment start times
    t_start = all_time[0]
    t_end = all_time[-1]
    t_grid = np.arange(t_start, t_end, t_step)

    for t0 in t_grid:
        # Select data in [t0, t0+dt]
        mask = (all_time >= t0) & (all_time <= t0 + dt)
        time = all_time[mask]
        flux = all_flux[mask]
        print(t0)

        # Skip if no data or gaps too large
        if len(time) == 0:
            continue
        if abs(time[-1] - t0 - dt) > expected_cadence / 86400:
            continue
        if abs(time[0] - t0) > expected_cadence / 86400:
            continue
        if np.any(np.diff(time) > expected_cadence / 86400):
            continue

        # Zero time for fitting
        t_zeroed = time - time[0]

        # Initial guesses for dt shift - tweak if needed
        dt_guesses = np.linspace(-0.5, 0.5, 10)

        t_est_list = []
        for dt_guess in dt_guesses:
            try:
                popt, _ = curve_fit(
                    lambda x, dt_shift: model_interp(x + dt_shift),
                    xdata=t_zeroed,
                    ydata=flux,
                    p0=[dt_guess],
                    xtol=1e-12,
                    maxfev=1000
                )
                # Convert estimated dt shift back to absolute time
                t_est_list.append(time[0] + popt[0])
            except Exception:
                pass

        if len(t_est_list) == 0:
            continue

        # Choose estimate closest to start time of segment
        t_est_arr = np.array(t_est_list)
        best_t_est = t_est_arr[np.argmin(np.abs(t_est_arr - time[0]))]

        true_time.append(time[0])
        est_time.append(best_t_est)

    true_time = np.array(true_time)
    est_time = np.array(est_time)

    # Find gaps to split into contiguous segments for plotting
    time_diff = np.diff(true_time)
    gap_mask = time_diff > 3000
    gap_indices = np.where(gap_mask)[0]
    segments = np.split(true_time, gap_indices + 1)

    # Compute t_shift for xlabel
    tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5) / 100) * 100)

    # Plotting
    margin = 0.5  # days
    fig_oc, axs = pt.subplots(
        1,
        len(segments),
        figsize=(8, 3),
        sharey=True,
        gridspec_kw={'wspace': 0, 'hspace': 0},
        width_ratios=[seg[-1] - seg[0] + margin * 2 for seg in segments]
    )

    if not isinstance(axs, np.ndarray):
        axs = [axs]

    for ii, ax in enumerate(axs):
        ax.scatter(
            true_time - tshift + OFFSET - 2400000.5,
            true_time - est_time,
            c='k',
            s=1
        )
        if ii == 0:
            ax.set_ylabel('O-C (days)', fontsize=11)
        ax.set_xlim(
            segments[ii][0] - tshift + OFFSET - 2400000.5 - margin,
            segments[ii][-1] - tshift + OFFSET - 2400000.5 + margin
        )

    fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
    pt.show()
    """

    """
    # Interpolate model profile
    model_phase = np.linspace(0, 1, len(model_profile), endpoint=False)
    interp_model = interp1d(model_phase, model_profile, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Prep light curve
    # === SETUP ===
    OFFSET = 2457000
    expected_cadence = 60  # seconds
    expected_cadence_days = expected_cadence / 86400
    t_step = 1    # days between segments
    dt = 1.0       # duration of each fit segment
    margin = 0.5     # for plot x-limits

    # === PREP MODEL ===
    model_phase = np.linspace(0, 1, len(model_profile), endpoint=False)
    interp_model = interp1d(model_phase, model_profile, kind='linear', bounds_error=False, fill_value='extrapolate')

    # === PREP LIGHT CURVE ===
    mask = np.isfinite(lc.flux.value.unmasked)
    all_flux = lc.flux.value[mask]
    all_time = lc.time.value[mask]

    # === FIND VALID SEGMENTS ===
    def find_valid_segments(all_time, all_flux, dt=1.0, t_step=0.1, expected_cadence_days=1800/86400):
        segments = []
        i = 0
        N = len(all_time)
        while i < N:
            print(i)
            t0 = all_time[i]
            t1 = t0 + dt
            j = np.searchsorted(all_time, t1)

            time_chunk = all_time[i:j]
            flux_chunk = all_flux[i:j]

            if len(time_chunk) >= 10:
                diffs = np.diff(time_chunk)
                if np.all(diffs < 1.5 * expected_cadence_days):  # allow small gaps
                    segments.append((time_chunk, flux_chunk))

            i = np.searchsorted(all_time, t0 + t_step)

        return segments

    valid_segments = find_valid_segments(all_time, all_flux, dt, t_step, expected_cadence_days)
    print("hellow")
    # === O–C CALCULATION ===
    period_grid = np.linspace(0.05, 0.2, 30)      # trial periods
    shift_grid = np.linspace(-0.2, 0.2, 60)       # trial phase shifts
    shift_grid_2D, period_grid_2D = np.meshgrid(shift_grid, period_grid)
    shift_flat = shift_grid_2D.ravel()
    period_flat = period_grid_2D.ravel()

    true_time = []
    est_time = []

    for time, flux in valid_segments:
        t0 = time[0]
        t_zeroed = time - t0
        N = len(t_zeroed)

        # Broadcast grid search
        time_matrix = t_zeroed[:, None] + shift_flat[None, :]
        phase_matrix = (time_matrix / period_flat[None, :]) % 1
        model_matrix = interp_model(phase_matrix)
        residuals = flux[:, None] - model_matrix
        mse = np.mean(residuals ** 2, axis=0)
        best_idx = np.argmin(mse)

        best_shift = shift_flat[best_idx]
        best_period = period_flat[best_idx]
        print(time)
        true_time.append(t0)
        est_time.append(t0 + best_shift)

    true_time = np.array(true_time)
    est_time = np.array(est_time)

    # === SEGMENT SPLITTING FOR PLOTTING ===
    time_diff = np.diff(true_time)
    gap_indices = np.where(time_diff > 3000 / 86400)[0]
    segments = np.split(true_time, gap_indices + 1)

    tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5) / 100) * 100)

    # === PLOT O–C DIAGRAM ===
    fig_oc, axs = pt.subplots(1, len(segments), figsize=(8, 3), sharey=True,
                            gridspec_kw={'wspace': 0, 'hspace': 0},
                            width_ratios=[seg[-1] - seg[0] + margin * 2 for seg in segments])

    if not isinstance(axs, np.ndarray):
        axs = [axs]

    for ii, ax in enumerate(axs):
        seg_mask = (true_time >= segments[ii][0]) & (true_time <= segments[ii][-1])
        oc_vals = true_time[seg_mask] - est_time[seg_mask]
        ax.scatter(true_time[seg_mask] - tshift + OFFSET - 2400000.5, oc_vals, c='k', s=1)
        if ii == 0:
            ax.set_ylabel('O-C (days)', fontsize=11)
        ax.set_xlim(segments[ii][0] - tshift + OFFSET - 2400000.5 - margin,
                    segments[ii][-1] - tshift + OFFSET - 2400000.5 + margin)

    fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
    pt.show()


    """
    return tshift


def get_csv_epsilon_value(csv_file_path): 

    try:
        df = pd.read_csv(csv_file_path)
        
        #if 'TIC_ID' not in df.columns:
        #    raise ValueError("CSV does not contain a 'TIC_ID' column.")
        
        KIC_list = df['KIC'].dropna().astype(str).tolist()
        FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
        i = 0 
        while( i < len(KIC_list)):
            get_epsilon_value(KIC_list[i], FUNCTION_list[i])
            i += 1
        
    except Exception as e:
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


#print(percenterror(-216, -258 ))

#getPeriodogramData('')
#print(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\tic_ids.csv"))

#fig, a = pt.subplots(1)
#getLightCurveData("KIC 3733346")
#identifyPeaks('X Caeli')
#identifyPeaks('X Caeli')
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

#seriesofstarsTest_time_error(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\KeplerStarsOutput_2_timeerror.csv"))
#identifyPeaksOfLightcurves_manual('KIC 3123138', 0)
#guessLegacy('KIC 4048494',0) 
#print(getMeanSquaredResidual('KIC 7548479'))
#identifyPeaks('KIC 12602250')
#print(guess("X Caeli", 1))
#3429637
get_csv_epsilon_value(r"ResearchPython\KeplerStarsOutput_combined.csv")
#lc = lk.search_lightcurve("KIC 3429637").download_all().stitch().remove_outliers(sigma = 5.0)
#pt.plot(lc.time.value, lc.flux.value)
#getPeriodogramData('KIC 3429637')s
#pt.title("Light Curve for KIC 3429637")
#pt.xlabel("Time (Days)")
#pt.ylabel("Flux")
#seriesofstarsTest(load_tic_ids_from_csv(r"C:\Users\ahmed\research_delta\tic_ids.csv"))
pt.show()

"""
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
