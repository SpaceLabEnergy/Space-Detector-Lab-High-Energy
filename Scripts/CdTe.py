
import gamma_ray_functions as gr_func

def main():
    # Energies of sources
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Am_energies = [59.5409]
    Cs_energies = [31.85, 661.657]
    
    # CdTe files
    CdTe_bkg_file = "CdTe_bck_10min.mca"
    CdTe_Ba_file = "CdTe_Ba_direct.mca"
    CdTe_Am_file = "CdTe_Am_direct.mca"
    CdTe_Cs_file = "CdTe_Cs_direct.mca"
    CdTe_Am_15_deg = "CdTe_Am_15.mca"
    CdTe_Am_30_deg = "CdTe_Am_30.mca"
    CdTe_Am_45_deg = "CdTe_Am_45.mca"
    CdTe_Am_60_deg = "CdTe_Am_60.mca"
    CdTe_Am_90_deg = "CdTe_Am_90.mca"
    
    # convert counts
    CdTe_Ba_channels, CdTe_Ba_count_rate, CdTe_Ba_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Ba_file, "CdTe")
    CdTe_Am_channels, CdTe_Am_count_rate, CdTe_Am_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_file, "CdTe")
    CdTe_Cs_channels, CdTe_Cs_count_rate, CdTe_Cs_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Cs_file, "CdTe")
    CdTe_Am15_channels, CdTe_Am15_count_rate, CdTe_Am15_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_15_deg, "CdTe")
    CdTe_Am30_channels, CdTe_Am30_count_rate, CdTe_Am30_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_30_deg, "CdTe")
    CdTe_Am45_channels, CdTe_Am45_count_rate, CdTe_Am45_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_45_deg, "CdTe")
    CdTe_Am60_channels, CdTe_Am60_count_rate, CdTe_Am60_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_60_deg, "CdTe")
    CdTe_Am90_channels, CdTe_Am90_count_rate, CdTe_Am90_count_rate_minus_bkg = gr_func.convert_counts(CdTe_bkg_file, CdTe_Am_90_deg, "CdTe")
    
    # fitting peaks
    CdTe_Ba_result, CdTe_Ba_centroids, CdTe_Ba_centroid_errors, CdTe_Ba_fwhm, CdTe_Ba_fwhm_errors = gr_func.fitting_peaks(CdTe_Ba_channels, CdTe_Ba_count_rate, 2, 2, (0.04, 0.25), None, None, "$Ba_{133}$", "CdTe", 1000)
    CdTe_Am_result, CdTe_Am_centroids, CdTe_Am_centroid_errors, CdTe_Am_fwhm, CdTe_Am_fwhm_errors = gr_func.fitting_peaks(CdTe_Am_channels, CdTe_Am_count_rate, 2, 2, 1.0, None, None, "$Am_{241}$", "CdTe", 1000)
    CdTe_Cs_result, CdTe_Cs_centroids, CdTe_Cs_centroid_errors, CdTe_Cs_fwhm, CdTe_Cs_fwhm_errors = gr_func.fitting_peaks(CdTe_Cs_channels, CdTe_Cs_count_rate, 5, 5, 0.06, None, None, "$Cs_{137}$", "CdTe", 600)
    CdTe_Am15_result, CdTe_Am15_centroids, CdTe_Am15_centroid_errors, CdTe_Am15_fwhm, CdTe_Am15_fwhm_errors = gr_func.fitting_peaks(CdTe_Am15_channels, CdTe_Am15_count_rate, 2, 2, 1.0, None, None, "$Am_{241}$ 15 deg", "CdTe", 1000)
    CdTe_Am30_result, CdTe_Am30_centroids, CdTe_Am30_centroid_errors, CdTe_Am30_fwhm, CdTe_Am30_fwhm_errors = gr_func.fitting_peaks(CdTe_Am30_channels, CdTe_Am30_count_rate, 2, 2, 1.0, None, None, "$Am_{241}$ 30 deg", "CdTe", 1000)
    CdTe_Am45_result, CdTe_Am45_centroids, CdTe_Am45_centroid_errors, CdTe_Am45_fwhm, CdTe_Am45_fwhm_errors = gr_func.fitting_peaks(CdTe_Am45_channels, CdTe_Am45_count_rate, 2, 2, 1.0, None, None, "$Am_{241}$ 45 deg", "CdTe", 1000)
    CdTe_Am60_result, CdTe_Am60_centroids, CdTe_Am60_centroid_errors, CdTe_Am60_fwhm, CdTe_Am60_fwhm_errors = gr_func.fitting_peaks(CdTe_Am60_channels, CdTe_Am60_count_rate, 2, 2, 0.8, None, None, "$Am_{241}$ 60 deg", "CdTe", 1000)
    CdTe_Am90_result, CdTe_Am90_centroids, CdTe_Am90_centroid_errors, CdTe_Am90_fwhm, CdTe_Am90_fwhm_errors = gr_func.fitting_peaks(CdTe_Am90_channels, CdTe_Am90_count_rate, 2, 2, 0.2, None, None, "$Am_{241}$ 90 deg", "CdTe", 1000)
    
    # define CdTe channels and energies and uncertainties
    CdTe_channels = CdTe_Ba_centroids + CdTe_Am_centroids + CdTe_Cs_centroids
    CdTe_channel_errors = CdTe_Ba_centroid_errors + CdTe_Am_centroid_errors + CdTe_Cs_centroid_errors
    CdTe_energies = Ba_energies + Am_energies + [Cs_energies[0]]
    
    # calibration curve
    CdTe_slope, CdTe_slope_err, CdTe_intercept, CdTe_intercept_err = gr_func.linear_fit_odr(CdTe_energies, CdTe_channels, CdTe_channel_errors, "CdTe")
    
    # define CdTe fwhm and uncertainties
    CdTe_fwhm = CdTe_Ba_fwhm + CdTe_Am_fwhm + CdTe_Cs_fwhm
    CdTe_fwhm_err = CdTe_Ba_fwhm_errors + CdTe_Am_fwhm_errors + CdTe_Cs_fwhm_errors
    
    # energy resolution
    CdTe_energy_resolution = gr_func.energy_resolution(CdTe_channels, CdTe_channel_errors, CdTe_fwhm, CdTe_fwhm_err, CdTe_slope, CdTe_slope_err, CdTe_intercept, CdTe_intercept_err, "CdTe")
    
    # absolute and intrinsic efficiency against energy
    CdTe_abs_eff_energy, CdTe_intr_eff_energy = gr_func.plot_CdTe_efficiencies(CdTe_Am_result, CdTe_Ba_result, CdTe_Cs_result)
    
    CdTe_source_results_array = [CdTe_Am_result, CdTe_Am15_result, CdTe_Am30_result, CdTe_Am45_result, CdTe_Am60_result, CdTe_Am90_result]
    angle_array = [0, 15, 30, 45, 60, 90]
    
    # absolute and intrinsic efficiency against angle
    CdTe_abs_eff_angle, CdTe_intr_eff_angle = gr_func.angle_efficiencies("CdTe", CdTe_source_results_array, angle_array)
    
if __name__ == "__main__":
    # there are no arguments to pass
    main()
