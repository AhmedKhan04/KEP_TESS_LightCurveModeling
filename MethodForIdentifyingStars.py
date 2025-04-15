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


def getLightCurveData(nameOfStar):
    search_result = lk.search_lightcurve(nameOfStar, quarter=(6,7,8))
    lc = search_result.download()
    #Modeling(l    lc.plot()
    
def getPeriodogramData(nameOfStar): 
    x = lk.search_targetpixelfile(nameOfStar).download()
    y = x.to_lightcurve()
    z = y.to_periodogram()
    z.smooth(method='logmedian', filter_width=0.1).plot(linewidth=2,  color='red', label='Smoothed', scale='log')
    z.plot(scale = 'log')
    #return z

def compGetPeriodogramData(nameOfStar): 
    x = lk.search_targetpixelfile(nameOfStar).download().to_lightcurve()
    y = lk.search_lightcurve(nameOfStar, quarter=(6,7,8)).download_all().stitch().remove_outliers(sigma = 5.0)
    
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

def identifyPeaks(nameOfStar):
    pg, lightc = compGetPeriodogramData(nameOfStar)
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=[max_power * 0.1, max_power * 1.1])
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
      
        flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
        amplitude_guess = flux_range
        phase_guess = 0  
        frequency_guess = frequencyfitted[b].value
        offset_guess = np.mean(flux)
        
        ig = [amplitude_guess, phase_guess, frequency_guess, offset_guess]
        
        # Adding bounds: to force some values
        #bounds = ([0.55*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        bounds = ([(bounds1/100)*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])

        if len(time) == 0 or len(flux) == 0:
            raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=10000, method='trf')
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
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=10000, method='trf')
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
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=10000, method='trf')
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
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
    amplitude_guess = flux_range
    phase_guess = 0 
    offset_guess = np.mean(flux)
    print(f"need to iterate: +  {len(frequencyfitted)} + times")

    while b < len(frequencyfitted):

        #Foldedlc = lc.fold(period = (1 / frequencyfitted[b].value))
 
        frequency_guess = frequencyfitted[b].value
        ig = [amplitude_guess, phase_guess, frequency_guess, offset_guess]
        x, y =  getMeanSquaredResidual(a,search_result, frequencyfitted, powers)
        # Adding bounds: to force some values of amplitude
        #bounds = ([0.55*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])
        bounds = ([y*amplitude_guess, -2*np.pi, 0.9*frequency_guess, np.min(flux)], [amplitude_guess, 2*np.pi, 1.1*frequency_guess, np.max(flux)])


        if len(time) == 0 or len(flux) == 0:
              raise ValueError("After cleaning, the time or flux array is empty.")
        #ig = [(np.max(flux) - np.min(flux))/2, 0, frequencyfitted[b].value, np.mean(flux)]
        params, _ = curve_fit(sine_model, time, flux, p0=ig, bounds=bounds, maxfev=10000, method='trf')
        amplitude, phase, frequency, offset = params
        fit_c = sine_model(time, *params)
        #amplitude, phase, frequency, offset = params
        #residuals = flux - (fit_c)
        #meanSquare = np.sum((residuals)**2)/len(flux)
        c.append(fit_c)
        params_list.append((amplitude, phase, frequency, offset)) 
        b += 1
        #print(f"Reduced Chi-squared: {reduced_chi_squared:.3f}")
    #print(f"Reduced Chi-squared Average: {np.mean(c):.3f}")
    return params_list, lc 

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
            sine_print_terms.append(f"{amplitude:.4f} * sin(2π * {frequency:.4f} * t + {phase:.4f}) + {offset:.4f}")
            p += 1
        #addedTogether  = addedTogether/total_weight
        print(f"Composite Sine Function for {a}:")
        print("f(t) = " + " + ".join(sine_print_terms))
        print(total_weight)
        return addedTogether, lc     
    

def plotsidebysideactual(a):
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
    Composite_function = getCompositeSine2(nameofStar) #2
    search_result = lk.search_lightcurve(nameofStar,quarter=(6,7,8))
    lc = search_result.download_all().stitch()
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
    np.set_printoptions(threshold=np.inf)
    print(peaksofcomposite)
    print(matched_lightcurve_peaks)
    matched_lightcurve_peaks_snipped = matched_lightcurve_peaks
    peaksofcomposite_snipped = peaksofcomposite
    print(f"Snipped peaks lengths: {len(peaksofcomposite_snipped)}, {len(matched_lightcurve_peaks_snipped)}")
    pt.figure(figsize=(12, 8))
    pt.scatter(peaksofcomposite_snipped, np.ones(peaksofcomposite_snipped.shape[0]), color='blue', label='composite')
    residuals = np.abs(peaksofcomposite_snipped - matched_lightcurve_peaks_snipped)
    print(residuals.shape)
    #pt.scatter(matched_lightcurve_peaks, np.array(fluxvaluesLight[:min_length]), color='red', label='Light Curve')
    pt.scatter(matched_lightcurve_peaks_snipped, np.ones(matched_lightcurve_peaks_snipped.shape[0]), color='red', label='Light Curve')
    pt.plot(matched_lightcurve_peaks_snipped,residuals,"o-", color = "black")
    print(residuals)
    print(f" The average residual in days is {np.mean(residuals)}")
    #pt.plot(time, flux[:min_length2], 'o-', color='black', label='Light Curve')
    #pt.plot(time, Composite_function[:min_length2], 'o-', color='green', label='Curve Fit')
    pt.tight_layout()
    pt.draw()
    pt.show()
    print("Original peaksofcomposite:", peaksofcomposite[:10])  # first 10 values
    print("Snipped peaksofcomposite:", peaksofcomposite_snipped[:10])  # first 10 values after snipping
    print("Original matched_lightcurve_peaks:", peaksoflightcurve[:10])
    print("Snipped matched_lightcurve_peaks:", matched_lightcurve_peaks_snipped[:10])
    print("StartingTime (index to cut):", startingTime)
                





#print(percenterror(-216, -258 ))

#getPeriodogramData('')

#fig, a = pt.subplots(1)
#getLightCurveData("KIC 3733346")
#identifyPeaks('X Caeli')
#identifyPeaks('X Caeli')
a = ['KIC 12602250', 'KIC 9700322','KIC 4048494', 'KIC 6951642', 'KIC 8623953', 'KIC 8197761']
b = ['KIC 3429637', 'KIC 10451090', 'KIC 2987660']
c = ['KIC 12602250' , 'KIC 9700322' , 'KIC 8197761' , 'KIC 8623953' ,  'KIC 6382916' ,'KIC 3429637']
d = ['KIC 2987660' , 'KIC 10451090' , 'KIC 8197761' , 'KIC 8623953' ,'KIC 3429637']

#KIC 8197761!!!!!!!!' 'V593 Lyr',
#getChiSquaredReduced('BO Lyn')
#print(getChiSquared('KIC 8197761'))
#plotsidebyside2('KIC 4733344') 
plotsidebysideactual('KIC 3429637')
#guessLegacy('KIC 4048494',0) 
#print(getMeanSquaredResidual('KIC 7548479'))
#identifyPeaks('KIC 12602250')
#print(guess("X Caeli", 1))
#3429637
#getPeriodogramData('KIC 3429637')s
#pt.title("Light Curve for KIC 3429637")
#pt.xlabel("Time (Days)")
#pt.ylabel("Flux")
pt.show()


