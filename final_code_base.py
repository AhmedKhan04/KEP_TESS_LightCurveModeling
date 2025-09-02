# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:04:19 2024

@author: ahmed
"""
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

def compGetPeriodogramData(nameOfStar): 
    """
    Helper function to get light curve and periodogram data.

    Args:
        nameOfStar (str): KIC code to search for.

    Returns:
        tuple: (Periodogram, LightCurve)
    """

    x = lk.search_targetpixelfile(nameOfStar).download().to_lightcurve()
    lightcurve = lk.search_lightcurve(nameOfStar).download_all().stitch().remove_outliers(sigma = 5.0)
    #y = lk.search_lightcurve(nameOfStar, quarter=(6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    periodogram = x.to_periodogram()
    return periodogram, lightcurve

def sine_model(t, amplitude, phase, frequency, offset):
    """
    Helper function to generate sine curve.

    Args:
        t (np.array): time array to generate sine curve for 

        amplitude (float): amplitude of sine curve 

        frequency (float): frequency of sine curve 

        offset (float): vertical shift of sine curve 
        
    Returns:
        sine_curve (np.array): array of applied sine curve to the time array.
    """

    sine_curve = amplitude * np.sin(2 * np.pi * frequency * t + phase) + offset
    return sine_curve

def identifyPeaks_dynamic(nameOfStar, lowerscalar = 0.1):
    """
    Helper function to dynamically identify peak frequencies of FFT.

    Args:
        nameOfStar (str): KIC code to search for.

        lowerscalar (int): scalar to set as minimum threshold of scalar * max_power required to identify peak

    Returns:
        freqs (np.array): frequencies of the peak frequinces identified of the FFT 
        
        lightc (Lightcurve): light curve of the star

        powers (np.array): array of power values for peak frequencies
    """
    pg, lightc = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * lowerscalar, max_power * 1.1])
    x = pg.frequency[peaks]
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
    freqs = pg.frequency[filtered_peaks],
    powers = pg.power[filtered_peaks]
    return(freqs, lightc, powers)

def identifyPeaks_non_dynamic(nameOfStar):
    """
    Helper function to non-dynamically efficiently identify peak frequencies of FFT.

    Args:
        nameOfStar (str): KIC code to search for.

    Returns:
        powers (np.array): array of power values for peak frequencies
        
        ltcurves (Lightcurve): light curve of the star
    """
    pg, ltcurves = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * 0.1, max_power * 1.1])
    x = pg.frequency[peaks]
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
    powers = pg.power[filtered_peaks]
    return(powers, ltcurves)

def guess_swift(nameOfStar,bounds1,search_result, frequencyfitted):
    """
    Generates a set of sinusoidal components for the predictive model swiftly

    Args:
       nameOfStar (str): KIC code to search for.

       bounds1 (float): initial amplitude bounding

       search_result (astroquery): astroquery of the lightcurve of the star

       frequencyfitted (np.array): peak frequencies of the stars pulsation 
    
    Returns:
        list_comp (np.array): list of sinusodial models for components of the predictive model
    
    """

    lc = search_result
    b = 0 
    list_comp = []
    while b < len(frequencyfitted):
        
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
        list_comp.append(fit_c)
        b = b + 1

    return list_comp

def guess_deep(a, scalar, frequencyfitted, search_result, powers):
    """
    Generates a set of sinusoidal components for the predictive model deeply

    Args:
       nameOfStar (str): KIC code to search for.

       scalar (float): initial amplitude bounding

       frequencyfitted (np.array): peak frequencies of the stars pulsation 

       search_result (astroquery): astroquery of the lightcurve of the star

       powers (np.array): array of power values for peak frequencies

    Returns:
        params_list (np.array): list of sinusoidal parameters for components of the predictive model

        lc (Lightcurve): lightcurve of the star

        c (np.array): list of sinusodial models for components of the predictive model
    
    """
    lc = search_result
    b = 0 
    list_comp = []
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

        list_comp.append(fit_c)

        amplitude, phase, frequency, offset = params
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1
        print("passed")
    return params_list, lc, list_comp

def align_arrays(time, flux):
    """
    Helper function to align time and flux arrays

    Args:
       time (np.array): time values array 

       flux (np.array): flux values array
    
    Returns:
        time (np.array): time values array aligned

        flux (np.array): flux values array aligned
    
    
    """
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    return time, flux

def getMeanSquaredResidual(nameOfStar, search_result, frequency, powerofpeaks_arg):
        """
        Helper function to calculate Mean Square Residual swiftly (Legacy)

        Args:
            nameOfStar (str): KIC code to search for.

            search_result (astroquery): astroquery of the lightcurve of the star

            frequency (np.array): peak frequencies of the stars pulsation 

            search_result (astroquery): astroquery of the lightcurve of the star

            powerofpeaks_arg (np.array): array of power values for peak frequencies

        Returns:
            bestmeanSquare (float): MSE Value

            bestbound (Lightcurve): bounding parameter for best results

            
        """
        bestmeanSquare = 100000
        bestBound = 0
        lc = search_result
        for bounds1 in range(54,56): 
            listofsines = guess_swift(nameOfStar,bounds1, search_result, frequency)
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
                
                total_weight += weight
                addedTogether += (weight * sinInterpolated)
                p += 1
                addedTogether  = addedTogether/total_weight
                residuals = flux - (addedTogether)
                meanSquare = np.sum((residuals)**2)/len(flux)
                if(meanSquare < bestmeanSquare):
                    bestmeanSquare = meanSquare
                    bestBound = bounds1
        print(bestmeanSquare)
        return bestmeanSquare, bestBound/100

def getResiduals(fit, flux): 
    """
        Helper function to calculate Residuals swiftly (Legacy)

        Args:
            fit (np.array): predictive model

            flux (np.array): flux values of the light curve
        
        Returns:
            meanSquare (float): MSE value 
    
    """
    residuals = flux - fit
    meanSquare = np.sum((residuals)**2)/len(flux)
    return meanSquare

def getCompositeSine2_swift(nameOfStar):
        """
        Swiftly generates the predictive model for star

        Args:
        
            nameOfStar (str): KIC code to search for.

        Returns:

            addedTogether (np.array): Predictive model for star
            
            lc (lightcurve): light curve object of star
        
        
        """
        powerOfPeaks, _ = identifyPeaks_non_dynamic(nameOfStar)
        print(len(powerOfPeaks))
        powerOfPeaks = powerOfPeaks.value
        listofsines, lc = guess_baseline(nameOfStar)
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

            sine_print_terms.append(f"{amplitude:.4f} * sin(2π * {frequency:.4f} * t + {phase:.4f}) + {offset:.4f}")
            p += 1
        print(f"Composite Sine Function for {nameOfStar}:")
        print("f(t) = " + " + ".join(sine_print_terms))
        print(total_weight)
        return addedTogether, lc     

def getCompositeSine2_deep(nameOfStar):
        """
        Deeply generates the predictive model for star

        Args:
        
            nameOfStar (str): KIC code to search for.

        Returns:

            addedTogether (np.array): Predictive model for star
            
            lc (lightcurve): light curve object of star

            composite_string: string of full predictive model as superpositions of sinusoidal functions
        
        """
        powerOfPeaks, _ = identifyPeaks_non_dynamic(nameOfStar)
        if(powerOfPeaks == -1):
            return [-10], 0
        print(len(powerOfPeaks))
        powerOfPeaks = powerOfPeaks.value
        frequencyfitted2, search_result2, powers2 = identifyPeaks_dynamic(nameOfStar)
        amplitude_scale = 0.5
        listofsines, lc, _ = guess_deep(nameOfStar, amplitude_scale, frequencyfitted2, search_result2, powers2)
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
            
            lower, _,fits_low = guess_deep(nameOfStar, low_amplitude_scale, frequencyfitted2, search_result2, powers2)
            upper, _, fits_high = guess_deep(nameOfStar, high_amplitude_scale, frequencyfitted2, search_result2, powers2)
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

            sine_print_terms.append(f"{amplitude:.4f} * sin(2π * {frequency:.4f} * t + {phase:.4f}) + {offset:.4f}")
            count += 1
        composite_string  = "f(t) = " + " + ".join(sine_print_terms)
        print(f"Composite Sine Function for {nameOfStar}:")
        print(composite_string)
        print(total_weight)
        print(listofindexs)
        return newaddedtogether, lc, composite_string    

def plotsidebyside_swift(nameOfStar):
    """
    Swiftly generates model and plots light curve, predictive model and residuals
    
    Args:
        nameOfStar (str): KIC code to search for.
    
    """
    function, lc = getCompositeSine2_swift(nameOfStar)
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
    print(residuals)
    pt.plot(time, residuals, 'o-', color='blue', label='O-C (Observed - Calculated)')
    pt.plot(time, flux, 'o-', color='red', label='Light Curve')
    pt.plot(time, function, 'o-', color='green', label='Curve Fit')
    #pt.plot(time, a, 'o-', color = 'blue')
    pt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    pt.title("O-C Diagram " + str(nameOfStar))
    pt.xlabel("Time (Days)")
    pt.ylabel("O-C (Flux Difference)")
    pt.legend()
    pt.grid()
    pt.tight_layout()
    pt.show()

def plotsidebyside_deep(nameOfStar):
    """
    Deeply generates model and plots light curve, predictive model and residuals
    
    Args:
        nameOfStar (str): KIC code to search for.
    
    """
    function, lc, _ = getCompositeSine2_deep(nameOfStar)
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
    print(residuals)
    pt.plot(time, residuals, 'o-', color='blue', label=' Residuals (Observed - Calculated)')
    pt.plot(time, flux, 'o-', color='red', label='Light Curve')
    pt.plot(time, function, 'o-', color='green', label='Curve Fit')
    #pt.plot(time, a, 'o-', color = 'blue')
    pt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    pt.title("Diagram For " + str(nameOfStar))
    pt.xlabel("Time -2454833 [BKJD Days]")
    pt.ylabel("Normalized Flux")
    pt.legend()
    pt.grid()
    pt.tight_layout()
    pt.show()

def interpolate(time, flux, target_time):
    """
    Helper function to interpolate values 

    Args:
        time (np.array): time array of light curve
        flux (np.array): flux array of light curve
        target_time (float): time point area desired
    
    Returns:
        interpolated_flux (np.array): interpolated flux array
    
    """
    ip = interp1d(time, flux, kind='nearest', bounds_error=False, fill_value='extrapolate')
    interpolated_flux = ip(target_time)
    return interpolated_flux 

def get_epsilon_value(star_name, sine_string):
    """
    For a given star and predictive model, generates epsilon values

    Args:
        star_name (str): name of star KIC

        sine_string (str): predictive model sine string
    
    Returns:
        epsilon_values (np.array): array containing all epsilon values
        standard_deviation (float): standard deviation of the residuals of epsilon values from the linear trend line
        slope_fit (float): slope of linear trendline of epsilon values
        dsct_per (float): period of the most dominant pulsation mode within the predictive model
    
    """

    search_result = lk.search_lightcurve(f"KIC {star_name}")
    lc = search_result.download_all().stitch().remove_outliers(sigma = 5.0)
    t = lc.time.value
    pattern = r'([+-]?\d*\.?\d+)\s*\*\s*sin\s*\(\s*2π\s*\*\s*([+-]?\d*\.?\d+)'  
    matches = re.findall(pattern, sine_string)
    if not matches:
        print(f"{star_name} did not work")
        return [-1], -1, -1, -1
    amp_freq_pairs = [(abs(float(amp)), float(freq)) for amp, freq in matches]
    max_amp, max_freq = max(amp_freq_pairs, key=lambda x: x[0])
    dsct_per =  1.0 / max_freq




    sine_string = sine_string.replace('sin', 'np.sin')
    sine_string = sine_string.replace('2π', '2 * np.pi ')
    #sine_string = sine_string.replace('t', ' * t ')
    sine_string = sine_string.replace("f(t) = ", "")


    OFFSET = 0
    expected_cadence = 1800  # seconds


    def create_model_function(sine_string):
        """Create a callable function from the sine string"""
        def model(t, dt, *params):
            shifted_t = t + dt + (OFFSET)

            print(sine_string)
            return eval(sine_string.replace('t', 'shifted_t'))
        return model

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
        
      
        t0_list = time[0] + np.linspace(-0.45,0.45,4) * dsct_per
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
   
    # plotting
    tshift = int(np.floor((true_time[0] + OFFSET - 2400000.5)/100)*100)
    margin = 0.5  # days

    """
    fig_oc, axs = pt.subplots(1, len(segments), figsize=(8,3), sharey=True, 
                            gridspec_kw={'wspace': 0, 'hspace': 0},
                            width_ratios=[seg[-1]-seg[0] + margin*2 for seg in segments])

    if not isinstance(axs, np.ndarray):
        axs = [axs]

    for ii, ax in enumerate(axs):
        ax.scatter(true_time-tshift + OFFSET - 2400000.5, true_time-est_time, c='k', s=0.5, label = "Phase Shift (BKJD Days)")
        if ii == 0:
            ax.set_ylabel(r'$\epsilon_{d}$ (Days)')
          
            pt.title(f"Epsilon Values for KIC {star_name}")
        ax.set_xlim(segments[ii][0]-tshift + OFFSET - 2400000.5 - margin, 
                    segments[ii][-1]-tshift + OFFSET - 2400000.5 + margin)
    
    positive_data = np.abs(true_time-est_time)
    m, b = np.polyfit(true_time-tshift + OFFSET - 2400000.5, positive_data, 1)
    x = true_time-tshift + OFFSET - 2400000.5
    pt.plot(x, x*m + b, color = 'r', label = "Linear Drift", linestyle = '--')
    fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)
    """
    data = true_time-est_time
    m, b = np.polyfit(true_time-tshift + OFFSET - 2400000.5, data, 1)
    x = true_time-tshift + OFFSET - 2400000.5

    #pt.plot(x, x*m + b, color = 'r', label = "Linear Drift", linestyle = '--')
    #fig_oc.supxlabel(f'Time (MJD) + {tshift}', fontsize=11, y=-0.05)



    regression = x*m + b 
    residuals = data - regression 
    sig = np.std(residuals)
    #print(f" normalized eps {np.mean(true_time-est_time)/dsct_per}")
    return true_time-est_time, sig, m, dsct_per

def guess_baseline(nameOfStar):
    """
    Given a star, generates a quick baseline predictive model guess

    Args:
        nameOfStar: KIC identifier for the star
    
    Returns:
        params_list (np.array): list of sinusoidal parameters for components of the predictive model

        lc (Lightcurve): lightcurve of the star
    
    """
    frequencyfitted, search_result, powers = identifyPeaks_dynamic(nameOfStar)
    lc = search_result
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
        c.append(fit_c)
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1

    return params_list, lc 

def getCompositeSine(nameOfStar):
        """
        Helper function to swiftly take sinusodial models and superimpose and scale them 

        Args:
            nameOfStar: KIC identifier of the star

        Return:
            addedTogether: predictive model of the star's lightcurve
        
        
        """
        listofsines = guess_baseline(nameOfStar)
        addedTogether = 0
        search_result = lk.search_lightcurve(nameOfStar,quarter=(6,7,8))
        lc = search_result.download_all().stitch()
        time = lc.time.value
        flux = lc.flux.value
        time, flux = align_arrays(time,flux)
        powerOfPeaks = identifyPeaks_non_dynamic(nameOfStar).value


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
        return addedTogether

def get_csv_epsilon_value(csv_file_path): 
    """
    For a given csv of stars and predictive models, generates epsilon values and saves to csv

    Args:
        csv_file_path (str): file path to data containing name of stars (KIC) and predictive model strings
    
    """
    print("Running")
    try:
        df = pd.read_csv(csv_file_path)
        KIC_list = df['KIC'].dropna().astype(str).tolist()
        FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
        i = 0 
        master_list_eps = []
        while( i < len(KIC_list)):
            eps, sig, m, P_max = get_epsilon_value(KIC_list[i], FUNCTION_list[i])
            if(i>0):
                pt.close('all')
            
            print(KIC_list[i])
            print(f"average eps {np.average(eps)}")
            print(f"standard dev {sig}")
            print(f"standard dev normalized {np.abs(sig/P_max)}")
            print(f"slope {m}")
            print(f"slope normalized {m/P_max}")
            master_list_eps.append({"KIC": KIC_list[i], "average eps": np.average(eps), "slope": m, "slope/P_MAX": m/P_max,"standard dev": sig, "standard dev/P_MAX": np.abs(sig/P_max)})
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
    """
    Helper function to devide segments of the light curve for epsilon vale calculations 

    Args:
        all_time (np.array): time array of light curve
        all_flux (np.array): flux array of light curve
        dt (float): size of window
        t_step (float): step size of window starting and ending points
        expected_cadence_days: cadence expected 

    Returns:
        segments (np.array): segment based array of all chunks needed for epsilon calculations
    
    
    """
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
                
def seriesofstars_deep(listofstars):
    """
    Given a list of stars generate deep predictive models and export MSE to csv

    Args:
        listofstars (np.array): KIC ID list of stars
    
    """
    results = []
    try:
        for star in listofstars:
            print(f"KIC {star}")
            function, lc, composite_strings = getCompositeSine2_deep(f"KIC {star}")
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

def SpectralResiduals(nameOfStar, sine_string): 
    """
    Given a star and its predictive model, generates the R2_LSP value between it and its light curve

    Args:
        nameOfStar (str): name of star desired in KIC
        
        sine_string (str): string consisting of full predicitve model as superpositions of sinusodial functions
    
    Returns:
        spec_res (float): spectral residuals between light curve and predictive model
        R2 (float): Normalized R2_LSP value between light curve and predictive model
    
    """
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

    def spectral_goodness_of_fit(time, lc, model):
        """
        Computes the spectral residual and normalized R^2_LSP goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - time: numpy array, the time array of the signal
        - lc: lightcurve object
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_LSP: float
        """
      
        lc_mod = lk.LightCurve(time=time, flux=model)
        pg_obs = lc.to_periodogram()
        pg_mod = lc_mod.to_periodogram(frequency = pg_obs.frequency)
        S_f = pg_obs.power.value/np.max(pg_obs.power.value)
        M_f = pg_mod.power.value/np.max(pg_mod.power.value)
        signal  = lc.flux.value/np.percentile(lc.flux.value, 90)
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
        S_bar = np.mean(signal)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_LSP = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_LSP
    
    spec_res, R2 = spectral_goodness_of_fit(t, lc, model)
    return spec_res, R2
    
def SpectralResidualsCsvBased(csv_file_path): 
    """
    Given csv of names of stars and their corresponding predictive model, generates the R2_LSP value between them and their light curve

    Args:
        csv_file_path (str): path to csv file input
    
    """
    
    print("Running")
    try:
        df = pd.read_csv(csv_file_path)
        

        
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
            master_list_eps.append({"KIC": KIC_list[i], "Spectral Residuals": spectral_resid, "R2_LSP": R2})
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

def plotMap():
    """
    Helper function to plot all stars analyzed using right assension and declination on map
    
    """
    pt.style.use(['science', 'no-latex'])

    with open(r"map_file_vizier", 'r') as file:
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
 
    pt.show()

def cleaning_tess(csv_path):
    """
    Given a csv containing KIC, TIC and predictive models for a series of stars, cleans the corresponding TESS data 
    and calculates the R2_LSP values. Saves it to an CSV output.

    Args:
        csv_path (str): path to the csv file as input
    
    """
    print("Running")

    df = pd.read_csv(csv_path)
    

    TIC_list = df['TIC'].dropna().astype(str).tolist()
    KIC_list = df['KIC'].dropna().astype(str).tolist()
    FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
    i = 0 
    master_lists_tess_pop = []
 
    def spectral_goodness_of_fit(signal, model):
        """
        Computes the spectral residual and normalized R^2_LSP goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - time: numpy array, the time array of the signal
        - lc: lightcurve object
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_LSP: float
        """

        time = np.array(time)
        signal = np.array(signal)
        model = np.array(model)
        mask = np.isfinite(signal) & np.isfinite(model) & np.isfinite(time)
        time = time[mask]
        signal = signal[mask]

        model = model[mask]
        lc_mod = lk.LightCurve(time = time, flux = model)
        lc_obs = lk.LightCurve(time = time, flux = signal)
        pg_obs = lc_obs.to_periodogram()
        pg_mod = lc_mod.to_periodogram(frequency = pg_obs.frequency)
        obs_power = np.array(pg_obs.power.value)
        mod_power = np.array(pg_mod.power.value)
        obs_power = np.nan_to_num(obs_power, nan=0.0)
        mod_power = np.nan_to_num(mod_power, nan=0.0)
        obs_power[obs_power < 0] = 0
        mod_power[mod_power < 0] = 0
        S_f = obs_power / np.max(obs_power)
        M_f = mod_power / np.max(mod_power)
        
        def FilterPeaksfft(fft, scalar = 0.4):
            fft /= np.max(fft)
            peaks, _= find_peaks(fft, prominence=np.max(fft) * scalar)
            if(len(peaks) == 0):
                print("using max")
                peaks = [np.argmax(fft)]
                print(peaks)
            peak_amps = fft[peaks]
            filtered_fft = np.zeros(len(fft))

            for peak_index in peaks:
                if peak_index < 100:
                    continue
                filtered_fft[peak_index - 25: peak_index+25] = fft[peak_index-25:peak_index+25]
            return filtered_fft
        S_f = FilterPeaksfft(S_f)
        M_f = FilterPeaksfft(M_f, scalar = 0.15)
        signal = (master_flux/np.percentile(master_flux, 90))

        pt.plot(np.arange(len(M_f)), S_f, label = 'Signal_clean')
        pt.plot(np.arange(len(M_f)), M_f, label = 'Model_clean')
        pt.title(f"KIC {KIC_list[i]} | TIC {TIC_list[i]}")
        pt.legend()
        #pt.show()
        pt.savefig(fr"{'path_file'}\pic_TIC_{TIC_list[i]}.png")
        pt.close()
       
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
       
      
        S_bar = np.mean(signal)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_LSP = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_LSP
    
    while( i < len(KIC_list)):
        try:
            master_flux = []
            master_time = []
            sine_string = FUNCTION_list[i]
            pt.close('all')

            result = lk.search_tesscut(f"TIC{TIC_list[i]}")
            search_result = result.table['sequence_number']
            
            print(search_result)
            tpf_collection = result.download_all(cutout_size=50)
            
            for l in tpf_collection:
                s = unpopular.Source(l.path, remove_bad=True)
                s.set_aperture(rowlims=[25, 26], collims=[25, 26])
                s.add_cpm_model(exclusion_size=5, n=64, predictor_method="similar_brightness")
                s.set_regs([0.1])
                s.holdout_fit_predict(k=100);

                aperture_normalized_flux = s.get_aperture_lc(data_type="normalized_flux")
                aperture_cpm_prediction = s.get_aperture_lc(data_type="cpm_prediction", weighting=None)

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
            master_flux = 2 * (master_flux - np.min(master_flux)) / (np.max(master_flux) - np.min(master_flux)) - 1
            spec, R2  = spectral_goodness_of_fit(master_time, master_flux, model)
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2LSP": R2})
            pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")
            i+=1

        except Exception as the_exception:
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": 0,"R2LSP": -1})
            print('did not find')
            print(the_exception)
            df = pd.DataFrame(master_lists_tess_pop)
            df.to_csv(fr"{'file_path'}\results.csv")
            i+=1
            continue

    df = pd.DataFrame(master_lists_tess_pop)
    df.to_csv(fr"{'file_path'}\results.csv")

# only for plotting  
def cleaning_tess_plotting(csv_path):
    """
    Given a csv containing KIC, TIC and predictive models for a series of stars, cleans the corresponding TESS data 
    and plots the FFT cleaning

    Args:
        csv_path (str): path to the csv file as input
    
    """
    print("Running")
    df = pd.read_csv(csv_path)

    TIC_list = df['TIC'].dropna().astype(str).tolist()
    KIC_list = df['KIC'].dropna().astype(str).tolist()
    FUNCTION_list = df['Composite Function'].dropna().astype(str).tolist()
    i = 0 
    master_lists_tess_pop = []
    

    def spectral_goodness_of_fit(time, signal, model):
        """
        Computes the spectral residual and normalized R^2_LSP goodness-of-fit
        between the signal and the model.
        
        Parameters:
        - time: numpy array, the time array of the signal
        - lc: lightcurve object
        - model: numpy array, the modeled signal
        
        Returns:
        - spectral_residual: float
        - R2_LSP: float
        """

        time = np.array(time)
        signal = np.array(signal)
        model = np.array(model)
        print(len(signal))
        mask = np.isfinite(signal) & np.isfinite(model) & np.isfinite(time)
        time = time[mask]
        
        signal = signal[mask]
        print(len(signal))
        model = model[mask]
        lc_mod = lk.LightCurve(time = time, flux = model)
        lc_obs = lk.LightCurve(time = time, flux = signal)
        pg_obs = lc_obs.to_periodogram()
        pg_mod = lc_mod.to_periodogram(frequency = pg_obs.frequency)
        obs_power = np.array(pg_obs.power.value)
        mod_power = np.array(pg_mod.power.value)
        obs_power = np.nan_to_num(obs_power, nan=0.0)
        mod_power = np.nan_to_num(mod_power, nan=0.0)
        obs_power[obs_power < 0] = 0
        mod_power[mod_power < 0] = 0
        S_f = obs_power / np.max(obs_power)
        M_f = mod_power / np.max(mod_power)
        pt.plot(pg_obs.frequency.value,S_f, color = 'red', label = 'signal')
        pt.plot(pg_obs.frequency.value,M_f, color = 'orange', label = 'model')
        
        def FilterPeaksfft(fft, scalar = 0.4):
            fft /= np.max(fft)
            peaks, _= find_peaks(fft, prominence=np.max(fft) * scalar)
            if(len(peaks) == 0):
                peaks = [np.argmax(fft)]
            peak_amps = fft[peaks]
            filtered_fft = np.zeros(len(fft))
            for peak_index in peaks:
                if peak_index < 100:
                    continue
                filtered_fft[peak_index - 25: peak_index+25] = fft[peak_index-25:peak_index+25]
            return filtered_fft
        S_f = FilterPeaksfft(S_f)
        M_f = FilterPeaksfft(M_f, scalar = 0.2)
        M_f = FilterPeaksfft(M_f, scalar = 0.15)
        signal = (master_flux/np.percentile(master_flux, 90))

        pt.plot(pg_obs.frequency.value, S_f, label = 'Signal_clean')
        pt.plot(pg_obs.frequency.value, M_f, label = 'Model_clean')
        pt.title(f"KIC {KIC_list[i]} | TIC {TIC_list[i]}")
        pt.xlabel('Frequency (cycles/BKJD)')
        pt.ylabel('Normalized Power')
        pt.legend()
        pt.show()
        spectral_residual = np.sum(np.abs(S_f - M_f)**2)
        S_bar = np.mean(signal)
        normalization = np.sum(np.abs(S_f - S_bar)**2)
        R2_LSP = 1 - (spectral_residual / normalization)
        
        return spectral_residual, R2_LSP
    
    while( i < len(KIC_list)):
        try:
            master_flux = []
            master_time = []
            sine_string = FUNCTION_list[i]
            pt.close('all')

            result = lk.search_tesscut(f"TIC{TIC_list[i]}")
            search_result = result.table['sequence_number']
            tpf_collection = result.download_all(cutout_size=50)

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
            master_flux = 2 * (master_flux - np.min(master_flux)) / (np.max(master_flux) - np.min(master_flux)) - 1
            spec, R2  = spectral_goodness_of_fit(master_time, master_flux, model )
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2LSP": R2})
            spec, R2  = spectral_goodness_of_fit(master_flux, model, master_time )
            print(spec, R2)
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": spec,"R2FFT": R2})
         
            pt.title(f"TIC {TIC_list[i]} | KIC {KIC_list[i]}")

            i+=1

        except Exception as the_exception:
            master_lists_tess_pop.append({"KIC": KIC_list[i],"TIC": TIC_list[i],"spectral_res": 0,"R2FFT": -1})
            print('did not find')
            print(the_exception)
            i+=1
            continue


        




if __name__ == "__main__":
    pt.style.use(['science', 'no-latex'])
    pt.rcParams.update({'figure.dpi': '300'})

    """
    Enter Code Here
    
    
    """

    pt.show()

