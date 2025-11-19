
"""
gamma_ray_functions.py

A collection of functions to help with data analysis
in the high energy detectors assignment.

"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Model
from lmfit.models import GaussianModel, LinearModel
from scipy.signal import find_peaks
from scipy import odr


def parse_spe_files(filename):
    """
    A function to parse a .SPE file containing headers and data.

    Parameters:
    -----------
    filename : str
        Path to the SPE file to parse.

    Returns:
    --------
    headers : dict
        A dictionary mapping header keys (e.g. "$DATE:", "$MEAS_TIM:") to
        their associated multiline string values.
    channels : numpy.ndarray
        An array of integer channel indices (0, 1, 2, ..., N-1).
    counts : numpy.ndarray
        An array of floating-point count values extracted from the $DATA section.
    """
    IS_DATA = False            # Tracks whether we've reached the data section
    headers = {}               # Dictionary to store the headers
    data = []                  # List to store the data

    current_key = None         # The current header key
    current_value_lines = []   # Gathers lines belonging to the current header

    with open(filename) as file:
        for line in file:
            line = line.strip()

            # detect start of data section
            if line == "$DATA:":
                # save the last header value before switching
                if current_key is not None:
                    headers[current_key] = "\n".join(current_value_lines).strip()
                IS_DATA = True
                continue

            # handle header lines
            if not IS_DATA:
                if line.startswith("$"):
                    # save previous header (if exists)
                    if current_key is not None:
                        headers[current_key] = "\n".join(current_value_lines).strip()
                    
                    # start a new header
                    current_key = line.strip()
                    current_value_lines = []
                else:
                    current_value_lines.append(line)
            else:
                # inside data section
                if line.startswith("0 1023"):
                    # skip known data line not wanted
                    continue
                try:
                    # convert data to floats
                    counts = float(line.strip())
                    data.append(counts)
                except ValueError:
                    # silently skip non-data lines
                    continue

    # in case file doesn’t end with $DATA:
    if not IS_DATA and current_key is not None:
        headers[current_key] = "\n".join(current_value_lines).strip()

    # build chanel and count arrays
    channels = np.arange(len(data))
    counts = np.array(data)

    return headers, channels, counts

def parse_mca_files(filename):
    """
    A function to parse an MCA file containing headers and data blocks.

    Parameters:
    -----------
    filename : str
        Path to the MCA file to be parsed.

    Returns:
    --------
    headers : dict
        Dictionary of header key/value pairs extracted from the
        "<<PMCA SPECTRUM>>" block.
    channels : numpy.ndarray
        Integer array of channel indices (0, 1, ..., N-1).
    counts : numpy.ndarray
        Integer count values extracted from the "<<DATA>>" block.
    """
    # Data flags to track header and data sections
    IS_HEADER = False
    IS_DATA = False
    
    # Set up lists for the header keys and values and data
    header_key = []
    header_value = []
    data = []

    with open(filename) as file:
        for line in file:
            line = line.strip()

            # Detect the beginning of the header block
            if line == "<<PMCA SPECTRUM>>":
                IS_HEADER = True
                IS_DATA = False
                continue

            # Detect the beginning of the data block
            elif line == "<<DATA>>":
                IS_HEADER = False
                IS_DATA = True
                continue

            # Detect any lines that aren't within the header and data blocks
            elif line.startswith("<<END") or line.startswith("<<DP5 CONFIGURATION") or line.startswith("<<DPP STATUS"):
                IS_HEADER = False
                IS_DATA = False
                continue

            # handle header lines
            if IS_HEADER:
                if "-" in line:
                    # Split the header lines by dashes
                    header_split = line.split("-", 1)
                    header_key.append(header_split[0].strip())
                    header_value.append(header_split[1].strip())
            
            # Handle data lines
            elif IS_DATA:
                try:
                    counts = int(line.strip())
                    data.append(counts)
                except ValueError:
                    # silently skip non-data lines
                    continue

    # build channel and count arrays
    channels = np.arange(len(data))
    counts = np.array(data)

    # create a dictionary for the headers using the lists of keys and values
    headers = {k: v for k, v in zip(header_key, header_value)}
    
    return headers, channels, counts

def counts_per_sec(counts, exposure_time):
    """
    A function to convert raw counts to count rate.
    
    Parameters:
    -----------
    counts: numpy.array
        Array of counts to be converted
    exposure_time: float
        Exposure time over which counts were collected
        
    Returns:
    --------
    counts_per_sec: numpy.array
        Count rate
    """
    # calculate count rate by dividing counts by exposure time
    counts_per_sec = counts / exposure_time
    return counts_per_sec

def convert_counts(bkg_file, file, detector):
    """
    A function to convert raw counts to counts-per-second 
    and apply background subtraction.

    Parameters
    ----------
    bkg_file : str
        Path to the background spectrum file.
    file : str
        Path to the spectrum file.
    detector : str
        Name of the detector, "BGO", "NaI", or "CdTe". 

    Returns
    -------
    channels : numpy.ndarray
        Channel indices corresponding to the spectrum.
    counts_sec : numpy.ndarray
        Count rate for the measured spectrum.
    counts_minus_bkg : numpy.ndarray
        Background-subtracted count rate.
    """

    # select parser function based on detector
    if "BGO" in detector or "NaI" in detector:
        # Parse source spectrum using (SPE file)
        headers, channels, counts = parse_spe_files(file)
        # Parse background spectrum (SPE file)
        bkg_headers, bkg_channels, bkg_counts = parse_spe_files(bkg_file)
    
        # Extract background exposure time
        meas_bkg = bkg_headers["$MEAS_TIM:"]
        time_bkg = meas_bkg.split(" ")
        bkg_exp_time = int(time_bkg[0])
    
        # Normalize background
        bkg_counts_per_sec = counts_per_sec(np.array(bkg_counts), bkg_exp_time)
    
        # Extract source exposure time
        meas_time = headers["$MEAS_TIM:"]
        time = meas_time.split(" ")
        exposure_time = int(time[0])
        
    else:
        # Parse source spectrum using (MCA file)
        headers, channels, counts = parse_mca_files(file)
        # Parse background spectrum (MCA file)
        bkg_headers, bkg_channels, bkg_counts = parse_mca_files(bkg_file)
    
        time_bkg = bkg_headers["REAL_TIME"]
        bkg_exp_time = float(time_bkg)
    
        bkg_counts_per_sec = counts_per_sec(np.array(bkg_counts), bkg_exp_time)
    
        time_data = headers["REAL_TIME"]
        exposure_time = float(time_data)
    
    # Convert raw counts to count rate
    counts_sec = counts_per_sec(np.array(counts), exposure_time)
    # Subtract background
    counts_minus_bkg = (np.array(counts_sec) - np.array(bkg_counts_per_sec)).flatten()
    
    return channels, counts_sec, counts_minus_bkg

def fitting_peaks(channels, counts, lower, upper, height, width, distance, title, detector, channel_max):
    """
    A function to detect and fit Gaussian peaks in a spectrum.

    Parameters
    ----------
    channels : list
        Channel numbers from the spectrum data.
    counts : list
        Count rate per channel.
    lower : int
        Number of channels to include to the left of each detected peak for fitting.
    upper : int
        Number of channels to include to the right of each detected peak for fitting.
    height : float
        Minimum peak height required for peak detection.
    width : float or int
        Minimum peak width required for peak detection.
    distance : float or int
        Minimum separation (in channels) between peaks.
    title : str
        Title of the plot and identifier for labeling results.
    detector : str
        Name of the detector ("NaI", "BGO", "CdTe", etc.) used to determine expected energies.
    channel_max : int
        Maximum x-axis limit for the plot.

    Returns
    -------
    fit_result : list
        List of lmfit ModelResult objects for each successfully fitted peak.
    centroids : list of float
        Best-fit Gaussian peak centers (in channel units).
    centroid_error : list of float
        Standard error estimates for each peak center.
    fwhm : list of float
        Best-fit full-width-at-half-maximum (FWHM) values.
    fwhm_error : list of float
        Standard error estimates for each FWHM.
    """
    # find peaks using scipy find_peaks function
    peaks, _ = find_peaks(counts, height=height, width=width, distance=distance)
    
    # Gaussian model for peak fits
    model = GaussianModel()
    
    # Lists for fit results and extracted parameters
    gauss_peaks = []
    centroids = []
    centroid_error = []
    fwhm = []
    fwhm_error = []
    fit_result = []
    
    # set up the plot
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(channels, counts)
    
    # Set up a counter for labelling fitted peaks
    peaknum = 1
    
    # Create a for loop to detect peaks individually
    for peak in peaks:
        # Define a fitting window around the peak
        left = max(0, peak - lower)
        right = min(len(channels), peak + upper)
        channel_range = np.arange(left, right)
        # Skip peaks too close to boundaries
        if peak < lower or peak + upper > len(channels):
            continue 
        
        # Initial parameter guess for Gaussian fit
        params = model.guess(counts[channel_range], channels[channel_range])
        
        # Perform the fit
        result = model.fit(counts[channel_range], x=channels[channel_range], params=params)
        
        # add results to arrays
        gauss_peaks.append(result)
        centroids.append(result.params["center"].value)
        centroid_error.append(result.params["center"].stderr)
        fwhm.append(result.params["fwhm"].value)
        fwhm_error.append(result.params["fwhm"].stderr)
        fit_result.append(result)
        
        # Plot initial and best fits
        ax.plot(channels[channel_range], result.init_fit, "--", label=f"initial fit peak {peaknum}")
        ax.plot(channels[channel_range], result.best_fit, "-", label=f"best fit peak {peaknum}")
        
        peaknum += 1
        
    # finalize spectrum plot
    ax.set_title(f"{title} Fitted Peaks {detector} Detector")
    ax.set_xlabel("Channels")
    ax.set_ylabel("Counts/sec")
    ax.legend()
    ax.set_xlim(0, channel_max)
    plt.show()
    
    # define expected gamma energies based on isotope name
    if "Ba" in title:
        energies = [53.1622, 80.9979, 356.0129]
    elif "Am" in title:
        energies = [59.5409]
    elif "Cs" in title:
        if "BGO" in detector or "NaI" in detector:
            energies = [661.657]
        elif "CdTe" in detector:
            energies = [31.85]
    else:
        energies = [1173.228, 1332.492]
        
    # Add text labels to the plot for each peak
    for energy, peak in zip(energies, peaks):
        plt.text(peak + 5, counts[peak] + 0.45, f"{energy:.1f} keV", color="red")
        
    # Print results
    print(f"\nGaussian Fit Results for {title}")
    print(f"--------------------------------")
    for cent, cent_error in zip(centroids, centroid_error):
        print(f"Channels      = {cent:.3f} ± {cent_error:.3f}")
    for f, f_err in zip(fwhm, fwhm_error):
        print(f"FWHM          = {f:.3f} ± {f_err:.3f}")
        
    return fit_result, centroids, centroid_error, fwhm, fwhm_error

def linear_fit_odr(energies, channels, channel_err, name):
    """
    Performs a linear calibration (Energy vs Channel) using scipy.odr,
    which properly accounts for uncertainties in the x-values (channels).

    Model: E = m*x + c

    Parameters
    ----------
    energies : list
        Known energies for calibration peaks (y-values)
    channels : list
        Channel numbers from the spectrum data.
    channel_err : list
        Uncertainties of the channel numbers
    name : str
        Name for plot labeling

    Returns
    -------
    output : scipy.odr.Output
        ODR fit result object
    """
    # Convert inputs to numpy arrays
    energies = np.asarray(energies, dtype=float)
    channels = np.asarray(channels, dtype=float)
    channel_err = np.asarray(channel_err, dtype=float)

    # Define linear model function
    def linear_func(B, x):
        return B[0] * x + B[1]

    # Create model and data for ODR
    model = odr.Model(linear_func)
    data = odr.RealData(channels, energies, sx=channel_err, sy=None)

    # Set up ODR with initial guesses
    beta0 = [1.0, 0.0]  # initial guess for slope, intercept
    odr_instance = odr.ODR(data, model, beta0=beta0)

    # Run the fit
    output = odr_instance.run()

    # Extract fit results
    slope, intercept = output.beta
    slope_err, intercept_err = output.sd_beta

    # Generate fitted line for plotting
    xfit = np.linspace(min(channels)*0.95, max(channels)*1.05, 300)
    yfit = linear_func(output.beta, xfit)

    # Plot results
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(channels, energies, xerr=channel_err, fmt='o', capsize=3,
                label='Data with x-errors')
    ax.plot(xfit, yfit, 'r-', label='ODR linear fit')
    ax.set_title(f"{name} Energies vs Channels (ODR Fit)")
    ax.set_xlabel(f"{name} Channels")
    ax.set_ylabel(f"{name} Energies (keV)")
    ax.legend()
    plt.show()

    # Print results
    print(f"\nLinear ODR Fit Results for {name}")
    print(f"--------------------------------")
    print(f"Slope      = {slope:.3f} ± {slope_err:.3f}")
    print(f"Intercept  = {intercept:.3f} ± {intercept_err:.3f}")
    print(f"Residual Variance = {output.res_var:.3f}")

    return slope, slope_err, intercept, intercept_err

def energy_resolution(channels, channel_errors, fwhm, fwhm_err, slope, slope_err, intercept, intercept_err, detector):
    """
    A function to compute the energy resolution of a detector and fit 
    it to a resolution model including uncertainty propagation.

    Parameters
    ----------
    channels : list
        Channel numbers from the spectrum data.
    channel_errors : numpy.array
        Uncertainties in the channel numbers.
    fwhm : list
        Full width at half maximum (FWHM) values from Gaussian fits (in channels).
    fwhm_err : list
        Uncertainties in FWHM values.
    slope : float
        Slope of energy calibration.
    slope_err : float
        Uncertainty in the calibration slope.
    intercept : float
        Energy calibration intercept.
    intercept_err : float
        Uncertainty in the calibration intercept.
    detector : str
        Name of the detector for figure labeling and printout.

    Returns
    -------
    None
        Prints results and produces a plot of energy resolution vs. energy.
    """
    
    # convert to arrays
    channels = np.asarray(channels, dtype=float)
    channel_errors = np.asarray(channel_errors, dtype=float)
    fwhm = np.asarray(fwhm, dtype=float)
    fwhm_err = np.asarray(fwhm_err, dtype=float)
    
    # convert channels to energies using calibration linear fit model
    def E_of_C(C):
        return slope*C + intercept

    energies = np.asarray(E_of_C(channels))

    # uncertainty propagation:
    E_err = np.sqrt(
        (channels * slope_err)**2 +
        (slope * channel_errors)**2 +
        intercept_err**2
    )
    
    # Compute ΔE and its uncertainty
    # derivative dE/dC = m
    dEdC = slope

    # ΔE = slope * FWHM
    fwhm_E = np.asarray(np.abs(fwhm * dEdC))

    # ΔE error propagation
    DeltaE_err = np.sqrt(
        (fwhm * slope_err)**2 +
        (slope * fwhm_err)**2
    )
    
    # Energy resolution and its uncertainty
    R = fwhm_E / energies

    R_err = R * np.sqrt(
        (DeltaE_err / fwhm_E)**2 +
        (E_err / energies)**2
    )

    R2 = R**2
    R2_err = 2 * R * R_err
    
    # Define energy resolution model
    def resolution_model(E, alpha, beta, gamma):
        return alpha*E**(-2) + beta*E**(-1) + gamma

    model = Model(resolution_model)
    
    # Initial guess
    params = model.make_params(alpha=1e-3, beta=1e-3, gamma=1e-3)

    # weights are inverse variance
    weights = 1.0 / R2_err

    # Fit the model
    result = model.fit(R2, params, E=energies, weights=weights)
    
    # Plot results
    plt.figure(figsize=(7,5))
    E_plot = np.linspace(min(energies), max(energies), 500)
    R2_fit_plot = model.eval(result.params, E=E_plot)

    plt.errorbar(energies, R, yerr=R_err, fmt='o', label="Data", color='blue')
    plt.plot(E_plot, np.sqrt(R2_fit_plot), color='red', label="Fit")

    plt.xlabel("Energy (keV)")
    plt.ylabel("Energy Resolution (ΔE/E)")
    plt.title(f"{detector} Energy Resolution Fit (with uncertainties)")
    plt.grid(True)
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    
    # Print results
    print(f"\nEnergy Resolution Fit Results for {detector}")
    print(f"--------------------------------")
    for R, R_err in zip(R, R_err):
        print(f"Energy Resolution      = {R:.3f} ± {R_err:.3f}")
    
    return

def find_efficiencies(detector, source, source_results):
    """
    A function to find absolute and intrinsic efficiencies of detector.
    
    Parameters:
    -----------
    detector: string 
        "BGO", "NaI", or "CdTe"
    source: string
        "AM", "BA", "CS", "CO"
    source_results: list
        Results from Gaussian fit of peaks in each source spectrum
        
    Returns:
    --------
    abs_eff: list
        Absolute efficiencies of detector
    intr_eff: list
        Intrinsic efficiencies of detector
    """

    ## Defining Source Activity and Detector Features
    # Source Activity Dec 1st 1979 (from box info)
    Am241_SA = 11.16 #microCurie
    Ba133_SA = 11.42
    Cs137_SA = 11.16
    Co60_SA = 12.28

    # Source Half-Life (from Student Handbook)
    Am241_HL = 432.2
    Ba133_HL = 10.537
    Cs137_HL = 30.09
    Co60_HL = 5.2714

    # Detector Diameter (from Datasheets)
    BGO_r = 2*2.54 # cm circular
    NaI_r = 2*2.54 # cm circular BUT NOT FROM CORRECT BIT OF DATASHEET
    CdTe_r = 3/10 # cm square

    # Detector Surface Area (cm**2)
    BGO_area = np.pi*(BGO_r**2)
    NaI_area = np.pi*(NaI_r**2)
    CdTe_area = CdTe_r**2

    # Source Distance to Detector (from Eris's measurements)
    BGO_d = 15 #cm
    NaI_d = 16 #cm
    CdTe_d = 10 #cm

    # Years since December 1st 1979 till Oct 21st (BGO) Oct 28th (NaI) and Nov 4 (CdTe)
    from datetime import datetime

    # Google Gemini generated code I turned into a function
    # (was going to just google: how long since __ but this is nicer)
    def how_long_yrs(date0, date1):
        from datetime import datetime

        # Calculate the difference
        difference = date1 - date0

        # Extract the years
        years = difference.days / 365
        
        return years

    # When we ran each detector
    date_start = datetime(1979, 12, 1)
    date_BGO = datetime(2025, 10, 21)
    date_NaI = datetime(2025, 10, 28)
    date_CdTe = datetime(2025, 11, 4)

    years_BGO = how_long_yrs(date_start, date_BGO)
    years_NaI = how_long_yrs(date_start, date_NaI)
    years_CdTe = how_long_yrs(date_start, date_CdTe)

    # Peak Energies
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Am_energies = [59.5409]
    Cs_energies = [661.657]
    Co_energies = [1173.228, 1332.492]

    # Create function that calculates current activity
    def activity_now(init_activity, half_life, time):
        import numpy as np
        
        current_activity = init_activity*np.exp(-(np.log(2)*time)/half_life)
        
        return current_activity

    ## Calculating Efficiencies

    abs_eff = []
    intr_eff = []


    if "Am" in source:
        # Americium has 1 Peak
        # Peak Total Counts Detected (in cps)
        Am_counts = source_results[0].params['amplitude'].value

        if "BGO" in detector:
            # Current Activity
            Am_CA = activity_now(Am241_SA, Am241_HL, years_BGO)*37000

            # Geometric Factor
            G0 = BGO_area / (4*np.pi*BGO_d**2)

        if "NaI" in detector:
            # Current Activity
            Am_CA = activity_now(Am241_SA, Am241_HL, years_NaI)*37000

            # Geometric Factor
            G0 = NaI_area / (4*np.pi*NaI_d**2)

        if "CdTe" in detector:
            # Current Activity
            Am_CA = activity_now(Am241_SA, Am241_HL, years_CdTe)*37000

            # Geometric Factor
            G0 = CdTe_area / (4*np.pi*CdTe_d**2)

        # Emmission Fraction
        Am_em_frac = 0.3578

        # Absolute Efficiency
        Am_abs_eff = Am_counts / (Am_CA * Am_em_frac)
        abs_eff.append(Am_abs_eff)

        # Intrinsic Efficiency
        Am_intr_eff = Am_counts / (Am_CA * G0 * Am_em_frac)
        intr_eff.append(Am_intr_eff)

    elif "Ba" in source:
        # Barium has 3 Peaks
        Ba_em_fracs = [0.0214, 0.329, 0.6205]
        index = 0
        
        if "BGO" in detector:
            # Current Activity
            Ba_CA = activity_now(Ba133_SA, Ba133_HL, years_BGO)*37000
        
            # Geometric Factor
            G0 = BGO_area / (4*np.pi*BGO_d**2)

        elif "NaI" in detector:
            # Current Activity
            Ba_CA = activity_now(Ba133_SA, Ba133_HL, years_NaI)*37000
            # Geometric Factor
            G0 = NaI_area / (4*np.pi*NaI_d**2)

        elif "CdTe" in detector:
            # Current Activity
            Ba_CA = activity_now(Ba133_SA, Ba133_HL, years_CdTe)*37000
            # Geometric Factor
            G0 = CdTe_area / (4*np.pi*CdTe_d**2)

        for result in source_results:
            # Peak Total Counts Detected (in cps)
            Ba_counts = result.params['amplitude'].value
        
            # Emmission Fraction
            Ba_em_frac = Ba_em_fracs[index]
        
            # Absolute Efficiency
            Ba_abs_eff = Ba_counts / (Ba_CA * Ba_em_frac)
            abs_eff.append(Ba_abs_eff)
        
            # Intrinsic Efficiency
            Ba_intr_eff = Ba_counts / (Ba_CA * G0 * Ba_em_frac)
            intr_eff.append(Ba_intr_eff)
            index += 1

    elif "Cs" in source:
        # Cesium has 1 Peak
        
        if "BGO" in detector:
            # Current Activity
            Cs_CA = activity_now(Cs137_SA, Cs137_HL, years_BGO)*37000
            # Geometric Factor
            G0 = BGO_area / (4*np.pi*BGO_d**2)

        elif "NaI" in detector:
            # Current Activity
            Cs_CA = activity_now(Cs137_SA, Cs137_HL, years_NaI)*37000
            # Geometric Factor
            G0 = NaI_area / (4*np.pi*NaI_d**2)

        elif "CdTe" in detector:
            # Current Activity
            Cs_CA = activity_now(Cs137_SA, Cs137_HL, years_CdTe)*37000
            # Geometric Factor
            G0 = CdTe_area / (4*np.pi*CdTe_d**2)

        # Peak Total Counts Detected (in cps)
        Cs_counts = source_results[0].params['amplitude'].value

        # Emmission Fraction
        Cs_em_frac = 0.8499

        # Absolute Efficiency
        Cs_abs_eff = Cs_counts / (Cs_CA * Cs_em_frac)
        abs_eff.append(Cs_abs_eff)

        # Intrinsic Efficiency
        Cs_intr_eff = Cs_counts / (Cs_CA * G0 * Cs_em_frac)
        intr_eff.append(Cs_intr_eff)

    elif "Co" in source:
        # Cobalt has 2 Peaks
        Co_em_fracs = [0.9985, 0.999826]
        index = 0

        if "BGO" in detector:
            # Current Activity
            Co_CA = activity_now(Co60_SA, Co60_HL, years_BGO)*37000

            #Geometric Factor
            G0 = BGO_area / (4*np.pi*BGO_d**2)

        elif "NaI" in detector:
            # Current Activity
            Co_CA = activity_now(Co60_SA, Co60_HL, years_NaI)*37000
            # Geometric Factor
            G0 = NaI_area / (4*np.pi*NaI_d**2)

        elif "CdTe" in detector:
            # Current Activity
            Co_CA = activity_now(Co60_SA, Co60_HL, years_CdTe)*37000
            # Geometric Factor
            G0 = CdTe_area / (4*np.pi*CdTe_d**2)

        for result in source_results:
            # Peak Total Counts Detected (in cps)
            Co_counts = result.params['amplitude'].value
            
            # Emmission Fraction
            Co_em_frac = Co_em_fracs[index]
        
            # Absolute Efficiency
            Co_abs_eff = Co_counts / (Co_CA * Co_em_frac)
            abs_eff.append(Co_abs_eff)
        
            # Intrinsic Efficiency
            Co_intr_eff = Co_counts / (Co_CA * G0 * Co_em_frac)
            intr_eff.append(Co_intr_eff)
        
            index += 1


    return abs_eff, intr_eff


def plot_efficiencies(detector, Am_results, Ba_results, Co_results, Cs_results):
    """
    A function to plot the detector absolute and intrinsic efficiencies
    as a function of energy.
    
    Parameters:
    detector: string
        "BGO" or "NaI"
    results: numpy.array
        result array of elements: source order: Am, Ba, Co, Cs (alphabetical)
        
    Returns:
    --------
    yval_a: list
        Absolute efficiencies
    yval_i: list
        Intrinsic efficiencies
    """
    # Peak Energies
    Am_energies = [59.5409]
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Co_energies = [1173.228, 1332.492]
    Cs_energies = [661.657]

    if "BGO" in detector:
        # Find Efficiencies
        Am_abs_eff, Am_intr_eff = find_efficiencies("BGO", "Am", Am_results)
        Ba_abs_eff, Ba_intr_eff = find_efficiencies("BGO", "Ba", Ba_results)
        Co_abs_eff, Co_intr_eff = find_efficiencies("BGO", "Co", Co_results)
        Cs_abs_eff, Cs_intr_eff = find_efficiencies("BGO", "Cs", Cs_results)

        # Title
        label_a = 'BGO Absolute Efficiency vs Energy'
        label_i = 'BGO Intrinsic Efficiency vs Energy'

    elif "NaI" in detector:
        # Find Efficiencies
        Am_abs_eff, Am_intr_eff = find_efficiencies("NaI", "Am", Am_results)
        Ba_abs_eff, Ba_intr_eff = find_efficiencies("NaI", "Ba", Ba_results)
        Co_abs_eff, Co_intr_eff = find_efficiencies("NaI", "Co", Co_results)
        Cs_abs_eff, Cs_intr_eff = find_efficiencies("NaI", "Cs", Cs_results)

        # Title
        label_a = 'NaI Absolute Efficiency vs Energy'
        label_i = 'NaI Intrinsic Efficiency vs Energy'

    # Absolute vs Energies
    yval_a = Am_abs_eff + Ba_abs_eff + Co_abs_eff + Cs_abs_eff
    xval = Am_energies + Ba_energies + Co_energies + Cs_energies

    fig, ax = plt.subplots(figsize=(9,5))
    ax.scatter(xval, yval_a)
    ax.set_xlabel(f'{detector} Peak Energies (keV)')
    ax.set_ylabel(f'{detector} Absolute Efficiency')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(label_a)


    # Intrinsic vs Energies
    yval_i = Am_intr_eff + Ba_intr_eff + Co_intr_eff + Cs_intr_eff

    fig, ax = plt.subplots(figsize=(9,5))
    ax.scatter(xval, yval_i)
    ax.set_xlabel(f'{detector} Peak Energies (keV)')
    ax.set_ylabel(f'{detector} Intrinsic Efficiency')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(label_i)
    
    print(f"\nAbsolute and Intrinsic Efficiencies vs Energy for {detector}")
    print(f"--------------------------------")
    for a, i in zip(yval_a, yval_i):
        print(f"Absolute Efficiency   = {a:.3f}")
        print(f"Intrinsic Efficiency  = {i:.3f}")
    
    return yval_a, yval_i

def plot_CdTe_efficiencies(Am_results, Ba_results, Cs_results):
    """
    A function to plot absolute and intrinsic efficiencies for the CdTe
    detector as a function of energy.
    
    Parameters:
    -----------
    results: numpy.array
        result array of elements: source order: Am, Ba, Co, Cs (alphabetical)
        
    Returns:
    --------
    yval_a: list
        Absolute efficiencies
    yval_i: list
        Intrinsic efficiencies
    """
    # Peak Energies
    Am_energies = [59.5409]
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Cs_energies = [661.657]

    # Find Efficiencies
    Am_abs_eff, Am_intr_eff = find_efficiencies("CdTe", "Am", Am_results)
    Ba_abs_eff, Ba_intr_eff = find_efficiencies("CdTe", "Ba", Ba_results)
    Cs_abs_eff, Cs_intr_eff = find_efficiencies("CdTe", "Cs", Cs_results)

    # Different Energy for Cesium
    Cs_energies = [31.85]

    # Title
    label_a = 'CdTe Absolute Efficiency vs Energy'
    label_i = 'CdTe Intrinsic Efficiency vs Energy'

    # Absolute vs Energies
    yval_a = Am_abs_eff + Ba_abs_eff + Cs_abs_eff
    xval = Am_energies + Ba_energies + Cs_energies

    fig, ax = plt.subplots(figsize=(9,5))
    ax.scatter(xval, yval_a)
    ax.set_xlabel('CdTe Peak Energies (keV)')
    ax.set_ylabel('CdTe Absolute Efficiency')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(label_a)

    # Intrinsic vs Energies
    yval_i = Am_intr_eff + Ba_intr_eff + Cs_intr_eff

    fig, ax = plt.subplots(figsize=(9,5))
    ax.scatter(xval, yval_i)
    ax.set_xlabel(f'CdTe Peak Energies (keV)')
    ax.set_ylabel(f'CdTe Intrinsic Efficiency')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(label_i)
    
    # Print results
    print(f"\nAbsolute and Intrinsic Efficiencies vs Energy for CdTe")
    print(f"--------------------------------")
    for a, i in zip(yval_a, yval_i):
        print(f"Absolute Efficiency   = {a:.3f}")
        print(f"Intrinsic Efficiency  = {i:.3f}")
    
    return yval_a, yval_i

def angle_efficiencies(detector, source_results_array, angle_array):
    """
    A function to find the absolute and intrinsic efficiencies of the 
    detector and plot them as a function of angle to investigate the
    off-axis performance of the detector.
    
    Parameters:
    -----------
    detector: string
        "BGO", "NaI", or "CdTe"
    source_results_array: numpy.array
        Array of the Gaussian fit results of the peaks for each source.
    angle_array: numpy.array
        Array of the angles the detector was set at.
        
    Returns:
    --------
    abs_eff: list
        Absolute efficiencies
    intr_eff: list
        Intrinsic efficiencies
    """
    # Americium has 1 Peak
    # Americium Values
    Am241_SA = 11.16 #microCurie
    Am241_HL = 432.2

    # Detector Diameter (from Datasheets)
    BGO_r = 2*2.54 # cm circular
    NaI_r = 2*2.54 # cm circular BUT NOT FROM CORRECT BIT OF DATASHEET
    CdTe_r = 3/10 # cm square

    # Detector Surface Area (cm**2)
    BGO_area = np.pi*(BGO_r**2)
    NaI_area = np.pi*(NaI_r**2)
    CdTe_area = CdTe_r**2

    # Source Distance to Detector (from Eris's measurements)
    BGO_d = 15 #cm
    NaI_d = 16 #cm
    CdTe_d = 10 #cm

    # Years since December 1st 1979 till Oct 21st (BGO) Oct 28th (NaI) and Nov 4 (CdTe)
    from datetime import datetime

    # Google Gemini generated code I turned into a function
    # (was going to just google: how long since __ but this is nicer)
    def how_long_yrs(date0, date1):
        from datetime import datetime
        # Calculate the difference
        difference = date1 - date0
        # Extract the years
        years = difference.days / 365
        
        return years

    # When we ran each detector
    date_start = datetime(1979, 12, 1)
    date_BGO = datetime(2025, 10, 21)
    date_NaI = datetime(2025, 10, 28)
    date_CdTe = datetime(2025, 11, 4)

    years_BGO = how_long_yrs(date_start, date_BGO)
    years_NaI = how_long_yrs(date_start, date_NaI)
    years_CdTe = how_long_yrs(date_start, date_CdTe)

    # Create function that calculates current activity
    def activity_now(init_activity, half_life, time):
        import numpy as np
        
        current_activity = init_activity*np.exp(-(np.log(2)*time)/half_life)
        
        return current_activity

    if "BGO" in detector:
        # Current Activity
        Am_CA = activity_now(Am241_SA, Am241_HL, years_BGO)*37000
        
        # Geometric Factor
        G0 = BGO_area / (4*np.pi*BGO_d**2)

    if "NaI" in detector:
        # Current Activity
        Am_CA = activity_now(Am241_SA, Am241_HL, years_NaI)*37000
        
        # Geometric Factor
        G0 = NaI_area / (4*np.pi*NaI_d**2)
        
    if "CdTe" in detector:
        # Current Activity
        Am_CA = activity_now(Am241_SA, Am241_HL, years_CdTe)*37000
        
        # Geometric Factor
        G0 = CdTe_area / (4*np.pi*CdTe_d**2)
        
    # Index for results array
    i = 0

    # Initialize lists, include direct observation
    #abs_eff_dir, intr_eff_dir = find_efficiencies("BGO", "AM", direct_results)
    abs_eff = []
    intr_eff = []

    for angle in angle_array:
        # Peak Total Counts Detected (in cps)
        results = source_results_array[i]
        Am_counts = results[0].params['amplitude'].value
        
        # Geometric Factor with Angle
        if angle == 90:
            G0_ang = G0
        else:
            G0_ang = np.cos((angle*np.pi)/180)*G0

        # Emmission Fraction
        Am_em_frac = 0.3578

        # Absolute Efficiency
        Am_abs_eff = Am_counts / (Am_CA * Am_em_frac)
        abs_eff.append(Am_abs_eff)

        # Intrinsic Efficiency
        Am_intr_eff = Am_counts / (Am_CA * G0_ang * Am_em_frac)
        intr_eff.append(Am_intr_eff)

        i += 1

    # Plot Angle vs Efficiencies
    fig, ax = plt.subplots(figsize=(9,5))
    ax.plot(angle_array, abs_eff, '-o')
    ax.set_xticks(np.arange(0,100,10))
    ax.set_xlabel(f'{detector} Detector Angle off Normal (deg)')
    ax.set_ylabel(f'{detector} Absolute Efficiency')
    ax.set_title(f'{detector} Americium Absolute Efficiency vs Angle')

    fig, ax = plt.subplots(figsize=(9,5))
    ax.plot(angle_array, intr_eff, '-o')
    ax.set_xticks(np.arange(0,100,10))
    ax.set_xlabel(f'{detector} Detector Angle off Normal (deg)')
    ax.set_ylabel(f'{detector} Intrinsic Efficiency')
    ax.set_title(f'{detector} Americium Intrinsic Efficiency vs Angle')
    
    print(f"\nAbsolute and Intrinsic Efficiencies vs Angle for {detector}")
    print(f"--------------------------------")
    for a, i in zip(abs_eff, intr_eff):
        print(f"Absolute Efficiency   = {a:.3f}")
        print(f"Intrinsic Efficiency  = {i:.3f}")
    
    return abs_eff, intr_eff
