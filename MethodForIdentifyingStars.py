import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lightkurve import search_lightcurve

def multi_sine_model(t, *params):
    """
    Model for multiple sine waves with different frequencies
    params format: [amp1, phase1, freq1, amp2, phase2, freq2, ..., offset]
    """
    n_sines = (len(params) - 1) // 3
    result = 0
    for i in range(n_sines):
        amp = params[i*3]
        phase = params[i*3 + 1]
        freq = params[i*3 + 2]
        result += amp * np.sin(2 * np.pi * freq * t + phase)
    
    # Add the offset (last parameter)
    result += params[-1]
    return result

def fit_multi_frequency(star_name, frequencies):
    """
    Fit multiple frequencies simultaneously to a light curve
    
    Parameters:
    -----------
    star_name : str
        Name of the star (e.g., 'KIC 3429637')
    frequencies : list
        List of frequencies to fit
        
    Returns:
    --------
    params : numpy array
        Fitted parameters
    lc : LightCurve object
        Light curve data
    """
    # Get light curve data
    search_result = search_lightcurve(star_name, quarter=(6,7,8))
    lc = search_result.download_all().stitch().remove_outliers(sigma=5.0)
    
    # Clean data
    time = lc.time.value
    flux = lc.flux.value
    vi = np.isfinite(time) & np.isfinite(flux)
    time = time[vi]
    flux = flux[vi]
    
    # Initial parameter guesses
    flux_range = np.percentile(flux, 95) - np.percentile(flux, 5)
    initial_params = []
    
    for freq in frequencies:
        # For each frequency, add amplitude, phase, and frequency
        amp_guess = flux_range * 0.1  # Start with a smaller amplitude guess
        phase_guess = 0
        freq_guess = freq
        initial_params.extend([amp_guess, phase_guess, freq_guess])
    
    # Add offset as the last parameter
    offset_guess = np.mean(flux)
    initial_params.append(offset_guess)
    
    # Set bounds - more relaxed than in your original code
    lower_bounds = []
    upper_bounds = []
    
    for freq in frequencies:
        # For each frequency
        lower_bounds.extend([0, -2*np.pi, 0.9*freq])  # Allow amplitudes to go to 0
        upper_bounds.extend([flux_range, 2*np.pi, 1.1*freq])
    
    # Bounds for offset
    lower_bounds.append(np.min(flux) - 0.1*flux_range)
    upper_bounds.append(np.max(flux) + 0.1*flux_range)
    
    # Curve fitting with all parameters at once
    try:
        params, _ = curve_fit(
            multi_sine_model, 
            time, 
            flux, 
            p0=initial_params,
            bounds=(lower_bounds, upper_bounds),
            maxfev=10000, 
            method='trf'  # Try different methods if this doesn't work well
        )
        
        # Calculate the fitted curve
        fit_curve = multi_sine_model(time, *params)
        
        # Calculate residuals and MSE
        residuals = flux - fit_curve
        mse = np.sum(residuals**2) / len(flux)
        print(f"MSE: {mse}")
        
        # Print the equation of the composite sine function
        print_multi_sine_equation(params, frequencies)
        
        return params, lc, fit_curve
        
    except Exception as e:
        print(f"Error in curve fitting: {e}")
        return None, lc, None

def print_multi_sine_equation(params, frequencies):
    """Print the equation of the multi-sine function with fitted parameters"""
    n_sines = (len(params) - 1) // 3
    sine_terms = []
    
    for i in range(n_sines):
        amp = params[i*3]
        phase = params[i*3 + 1]
        freq = params[i*3 + 2]
        sine_terms.append(f"{amp:.6f} * sin(2π * {freq:.6f} * t + {phase:.6f})")
    
    offset = params[-1]
    equation = " + ".join(sine_terms) + f" + {offset:.6f}"
    print(f"Composite Sine Function:")
    print(f"f(t) = {equation}")

def plot_results(star_name, time, flux, fit_curve):
    """Plot the original light curve, fitted curve, and residuals"""
    residuals = flux - fit_curve
    
    plt.figure(figsize=(12, 8))
    
    # Plot light curve and fit
    plt.subplot(2, 1, 1)
    plt.plot(time, flux, 'o', color='red', alpha=0.5, markersize=2, label='Light Curve')
    plt.plot(time, fit_curve, '-', color='green', linewidth=1, label='Curve Fit')
    plt.title(f"Light Curve and Fit for {star_name}")
    plt.xlabel("Time (Days)")
    plt.ylabel("Flux")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot residuals
    plt.subplot(2, 1, 2)
    plt.plot(time, residuals, 'o', color='blue', alpha=0.5, markersize=2, label='O-C (Observed - Calculated)')
    plt.axhline(0, color='red', linestyle='--', linewidth=1, label='Zero Line')
    plt.title(f"O-C Diagram for {star_name}")
    plt.xlabel("Time (Days)")
    plt.ylabel("O-C (Flux Difference)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def identify_frequencies(star_name, min_freq=1.0, height_factor=0.1, max_peaks=5):
    """
    Identify dominant frequencies from periodogram
    
    Parameters:
    -----------
    star_name : str
        Name of the star
    min_freq : float
        Minimum frequency to consider
    height_factor : float
        Factor to determine peak height threshold
    max_peaks : int
        Maximum number of peaks to return
        
    Returns:
    --------
    list of frequencies
    """
    from scipy.signal import find_peaks
    import lightkurve as lk
    
    # Get periodogram
    tpf = lk.search_targetpixelfile(star_name).download()
    lc = tpf.to_lightcurve()
    pg = lc.to_periodogram()
    
    # Find peaks
    max_power = np.max(pg.power.value)
    peaks, _ = find_peaks(pg.power, height=max_power * height_factor)
    
    # Filter peaks
    frequencies = []
    for i in range(len(peaks)):
        freq_value = pg.frequency[peaks[i]].value
        
        # Only consider frequencies above minimum
        if freq_value >= min_freq:
            # Check if similar to any existing frequency
            is_unique = True
            for existing_freq in frequencies:
                if abs(freq_value - existing_freq) <= 0.3:  # Tolerance for similar frequencies
                    is_unique = False
                    break
            
            if is_unique:
                frequencies.append(freq_value)
        
        # Limit number of frequencies
        if len(frequencies) >= max_peaks:
            break
    
    print(f"Identified frequencies: {frequencies}")
    return frequencies

def analyze_star(star_name):
    """Complete analysis of a star"""
    # Identify frequencies
    frequencies = identify_frequencies(star_name)
    
    if not frequencies:
        print("No significant frequencies found")
        return
    
    # Fit multi-frequency model
    params, lc, fit_curve = fit_multi_frequency(star_name, frequencies)
    
    if params is not None:
        # Get time and flux for plotting
        time = lc.time.value
        flux = lc.flux.value
        vi = np.isfinite(time) & np.isfinite(flux)
        time = time[vi]
        flux = flux[vi]
        
        # Plot results
        plot_results(star_name, time, flux, fit_curve)
    else:
        print("Fitting failed")

# Example usage
if __name__ == "__main__":
    # Analyze a single star
    analyze_star('KIC 3429637')
    
    # Or analyze multiple stars
    # stars = ['KIC 12602250', 'KIC 9700322', 'KIC 8197761']
    # for star in stars:
    #     print(f"\nAnalyzing {star}...")
    #     analyze_star(star)