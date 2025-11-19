
import gamma_ray_functions as gr_func

def main():
    # Energies of sources
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Am_energies = [59.5409]
    Cs_energies = [31.85, 661.657]
    Co_energies = [1173.228, 1332.492]
    
    # BGO files
    BGO_bkg_file = "BGO_Bck_10min.Spe"
    BGO_Ba_file = "BGO_Ba_direct.Spe"
    BGO_Am_file = "BGO_Am_241_direct.Spe"
    BGO_Cs_file = "BGO_Cs137_direct.Spe"
    BGO_Co_file = "bgo-co60-000-20cm-u.Spe"
    BGO_Am_15_deg = "BGO_Am_241_15deg.Spe"
    BGO_Am_30_deg = "BGO_Am_241_30deg.Spe"
    BGO_Am_45_deg = "BGO_Am_241_45deg.Spe"
    BGO_Am_60_deg = "BGO_Am241_60deg.Spe"
    BGO_Am_90_deg = "BGO_Am_241_90deg.Spe"
    
    # convert counts
    BGO_Ba_channels, BGO_Ba_count_rate, BGO_Ba_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Ba_file, "BGO")
    BGO_Am_channels, BGO_Am_count_rate, BGO_Am_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_file, "BGO")
    BGO_Cs_channels, BGO_Cs_count_rate, BGO_Cs_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Cs_file, "BGO")
    BGO_Co_channels, BGO_Co_count_rate, BGO_Co_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Co_file, "BGO")
    BGO_Am15_channels, BGO_Am15_count_rate, BGO_Am15_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_15_deg, "BGO")
    BGO_Am30_channels, BGO_Am30_count_rate, BGO_Am30_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_30_deg, "BGO")
    BGO_Am45_channels, BGO_Am45_count_rate, BGO_Am45_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_45_deg, "BGO")
    BGO_Am60_channels, BGO_Am60_count_rate, BGO_Am60_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_60_deg, "BGO")
    BGO_Am90_channels, BGO_Am90_count_rate, BGO_Am90_count_rate_minus_bkg = gr_func.convert_counts(BGO_bkg_file, BGO_Am_90_deg, "BGO")
    
    # fitting peaks
    BGO_Ba_result, BGO_Ba_centroids, BGO_Ba_centroid_errors, BGO_Ba_fwhm, BGO_Ba_fwhm_errors = gr_func.fitting_peaks(BGO_Ba_channels, BGO_Ba_count_rate_minus_bkg, 10, 20, 0.5, 5, None, "$Ba_{133}$", "BGO", 400)
    BGO_Am_result, BGO_Am_centroids, BGO_Am_centroid_errors, BGO_Am_fwhm, BGO_Am_fwhm_errors = gr_func.fitting_peaks(BGO_Am_channels, BGO_Am_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$", "BGO", 400)
    BGO_Cs_result, BGO_Cs_centroids, BGO_Cs_centroid_errors, BGO_Cs_fwhm, BGO_Cs_fwhm_errors = gr_func.fitting_peaks(BGO_Cs_channels, BGO_Cs_count_rate_minus_bkg, 10, 20, 5.0, 20, None, "$Cs_{137}$", "BGO", 600)
    BGO_Co_result, BGO_Co_centroids, BGO_Co_centroid_errors, BGO_Co_fwhm, BGO_Co_fwhm_errors = gr_func.fitting_peaks(BGO_Co_channels, BGO_Co_count_rate_minus_bkg, 20, 20, (0.05, 0.07), 6, 5, "$Co_{60}$", "BGO", max(BGO_Co_channels))
    BGO_Am15_result, BGO_Am15_centroids, BGO_Am15_centroid_errors, BGO_Am15_fwhm, BGO_Am15_fwhm_errors = gr_func.fitting_peaks(BGO_Am15_channels, BGO_Am15_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$ 15 deg", "BGO", 400)
    BGO_Am30_result, BGO_Am30_centroids, BGO_Am30_centroid_errors, BGO_Am30_fwhm, BGO_Am30_fwhm_errors = gr_func.fitting_peaks(BGO_Am30_channels, BGO_Am30_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$ 30 deg", "BGO", 400)
    BGO_Am45_result, BGO_Am45_centroids, BGO_Am45_centroid_errors, BGO_Am45_fwhm, BGO_Am45_fwhm_errors = gr_func.fitting_peaks(BGO_Am45_channels, BGO_Am45_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$ 45 deg", "BGO", 400)
    BGO_Am60_result, BGO_Am60_centroids, BGO_Am60_centroid_errors, BGO_Am60_fwhm, BGO_Am60_fwhm_errors = gr_func.fitting_peaks(BGO_Am60_channels, BGO_Am60_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$ 60 deg", "BGO", 400)
    BGO_Am90_result, BGO_Am90_centroids, BGO_Am90_centroid_errors, BGO_Am90_fwhm, BGO_Am90_fwhm_errors = gr_func.fitting_peaks(BGO_Am90_channels, BGO_Am90_count_rate_minus_bkg, 10, 20, 20, None, None, "$Am_{241}$ 90 deg", "BGO", 400)
    
    # define BGO channels and energies and uncertainties
    BGO_channels = BGO_Ba_centroids + BGO_Am_centroids + BGO_Cs_centroids + BGO_Co_centroids
    BGO_channel_errors = BGO_Ba_centroid_errors + BGO_Am_centroid_errors + BGO_Cs_centroid_errors + BGO_Co_centroid_errors
    BGO_energies = Ba_energies + Am_energies + [Cs_energies[1]] + Co_energies
    
    # calibration curve
    BGO_slope, BGO_slope_err, BGO_intercept, BGO_intercept_err = gr_func.linear_fit_odr(BGO_energies, BGO_channels, BGO_channel_errors, "BGO")
    
    # define BGO fwhm and uncertainties
    BGO_fwhm = BGO_Ba_fwhm + BGO_Am_fwhm + BGO_Cs_fwhm + BGO_Co_fwhm
    BGO_fwhm_err = BGO_Ba_fwhm_errors + BGO_Am_fwhm_errors + BGO_Cs_fwhm_errors + BGO_Co_fwhm_errors
    
    # energy resolution
    BGO_energy_resolution = gr_func.energy_resolution(BGO_channels, BGO_channel_errors, BGO_fwhm, BGO_fwhm_err, BGO_slope, BGO_slope_err, BGO_intercept, BGO_intercept_err, "BGO")
    
    # absolute and intrinsic efficiency against energy
    BGO_abs_eff_energy, BGO_intr_eff_energy = gr_func.plot_efficiencies("BGO", BGO_Am_result, BGO_Ba_result, BGO_Co_result, BGO_Cs_result)
    
    BGO_source_results_array = [BGO_Am_result, BGO_Am15_result, BGO_Am30_result, BGO_Am45_result, BGO_Am60_result, BGO_Am90_result]
    angle_array = [0, 15, 30, 45, 60, 90]
    
    # absolute and intrinsic efficiency against angle
    BGO_abs_eff_angle, BGO_intr_eff_angle = gr_func.angle_efficiencies("BGO", BGO_source_results_array, angle_array)
    
if __name__ == "__main__":
    # there are no arguments to pass
    main()
