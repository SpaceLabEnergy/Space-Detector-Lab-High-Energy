
import gamma_ray_functions as gr_func

def main():
    # Energies of sources
    Ba_energies = [53.1622, 80.9979, 356.0129]
    Am_energies = [59.5409]
    Cs_energies = [31.85, 661.657]
    Co_energies = [1173.228, 1332.492]
    
    # NaI files
    NaI_bkg_file = "NaI_bck_10min.Spe"
    NaI_Ba_file = "NaI_Ba_direct.Spe"
    NaI_Am_file = "NaI_Am_direct.Spe"
    NaI_Cs_file = "NaI_Cs_direct.Spe"
    NaI_Co_file = "NaI-co60-000-10cm-u.Spe"
    NaI_Am_15_deg = "NaI_Am_15deg.Spe"
    NaI_Am_30_deg = "NaI_Am_30deg.Spe"
    NaI_Am_45_deg = "NaI_Am_45deg.Spe"
    NaI_Am_60_deg = "NaI_Am_60deg.Spe"
    NaI_Am_90_deg = "NaI_Am_90deg.Spe"
    
    # convert counts
    NaI_Ba_channels, NaI_Ba_count_rate, NaI_Ba_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Ba_file, "NaI")
    NaI_Am_channels, NaI_Am_count_rate, NaI_Am_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_file, "NaI")
    NaI_Cs_channels, NaI_Cs_count_rate, NaI_Cs_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Cs_file, "NaI")
    NaI_Co_channels, NaI_Co_count_rate, NaI_Co_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Co_file, "NaI")
    NaI_Am15_channels, NaI_Am15_count_rate, NaI_Am15_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_15_deg, "NaI")
    NaI_Am30_channels, NaI_Am30_count_rate, NaI_Am30_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_30_deg, "NaI")
    NaI_Am45_channels, NaI_Am45_count_rate, NaI_Am45_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_45_deg, "NaI")
    NaI_Am60_channels, NaI_Am60_count_rate, NaI_Am60_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_60_deg, "NaI")
    NaI_Am90_channels, NaI_Am90_count_rate, NaI_Am90_count_rate_minus_bkg = gr_func.convert_counts(NaI_bkg_file, NaI_Am_90_deg, "NaI")
    
    # fitting peaks
    NaI_Ba_result, NaI_Ba_centroids, NaI_Ba_centroid_errors, NaI_Ba_fwhm, NaI_Ba_fwhm_errors = gr_func.fitting_peaks(NaI_Ba_channels, NaI_Ba_count_rate_minus_bkg, 10, 20, (0.5, 6), 3, 10, "$Ba_{133}$", "NaI", 400)
    NaI_Am_result, NaI_Am_centroids, NaI_Am_centroid_errors, NaI_Am_fwhm, NaI_Am_fwhm_errors = gr_func.fitting_peaks(NaI_Am_channels, NaI_Am_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$", "NaI", 400)
    NaI_Cs_result, NaI_Cs_centroids, NaI_Cs_centroid_errors, NaI_Cs_fwhm, NaI_Cs_fwhm_errors = gr_func.fitting_peaks(NaI_Cs_channels, NaI_Cs_count_rate_minus_bkg, 10, 20, 4.5, 5, None, "$Cs_{137}$", "NaI", 600)
    NaI_Co_result, NaI_Co_centroids, NaI_Co_centroid_errors, NaI_Co_fwhm, NaI_Co_fwhm_errors = gr_func.fitting_peaks(NaI_Co_channels, NaI_Co_count_rate, 10, 20, (0.01, 0.25), 6.9, 10, "$Co_{60}$", "NaI", max(NaI_Co_channels))
    NaI_Am15_result, NaI_Am15_centroids, NaI_Am15_centroid_errors, NaI_Am15_fwhm, NaI_Am15_fwhm_errors = gr_func.fitting_peaks(NaI_Am15_channels, NaI_Am15_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$ 15 deg", "NaI", 400)
    NaI_Am30_result, NaI_Am30_centroids, NaI_Am30_centroid_errors, NaI_Am30_fwhm, NaI_Am30_fwhm_errors = gr_func.fitting_peaks(NaI_Am30_channels, NaI_Am30_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$ 30 deg", "NaI", 400)
    NaI_Am45_result, NaI_Am45_centroids, NaI_Am45_centroid_errors, NaI_Am45_fwhm, NaI_Am45_fwhm_errors = gr_func.fitting_peaks(NaI_Am45_channels, NaI_Am45_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$ 45 deg", "NaI", 400)
    NaI_Am60_result, NaI_Am60_centroids, NaI_Am60_centroid_errors, NaI_Am60_fwhm, NaI_Am60_fwhm_errors = gr_func.fitting_peaks(NaI_Am60_channels, NaI_Am60_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$ 60 deg", "NaI", 400)
    NaI_Am90_result, NaI_Am90_centroids, NaI_Am90_centroid_errors, NaI_Am90_fwhm, NaI_Am90_fwhm_errors = gr_func.fitting_peaks(NaI_Am90_channels, NaI_Am90_count_rate_minus_bkg, 10, 20, 0.5, 3, None, "$Am_{241}$ 90 deg", "NaI", 400)
    
    # define NaI channels and energies and uncertainties
    NaI_channels = NaI_Ba_centroids + NaI_Am_centroids + NaI_Cs_centroids + NaI_Co_centroids
    NaI_channel_errors = NaI_Ba_centroid_errors + NaI_Am_centroid_errors + NaI_Cs_centroid_errors + NaI_Co_centroid_errors
    NaI_energies = Ba_energies + Am_energies + [Cs_energies[1]] + Co_energies
    
    # calibration curve
    NaI_slope, NaI_slope_err, NaI_intercept, NaI_intercept_err = gr_func.linear_fit_odr(NaI_energies, NaI_channels, NaI_channel_errors, "NaI")
    
    # define NaI fwhm and uncertainties
    NaI_fwhm = NaI_Ba_fwhm + NaI_Am_fwhm + NaI_Cs_fwhm + NaI_Co_fwhm
    NaI_fwhm_err = NaI_Ba_fwhm_errors + NaI_Am_fwhm_errors + NaI_Cs_fwhm_errors + NaI_Co_fwhm_errors
    
    # energy resolution
    NaI_energy_resolution = gr_func.energy_resolution(NaI_channels, NaI_channel_errors, NaI_fwhm, NaI_fwhm_err, NaI_slope, NaI_slope_err, NaI_intercept, NaI_intercept_err, "NaI")
    
    # absolute and intrinsic efficiency against energy
    NaI_abs_eff_energy, NaI_intr_eff_energy = gr_func.plot_efficiencies("NaI", NaI_Am_result, NaI_Ba_result, NaI_Co_result, NaI_Cs_result)
    
    NaI_source_results_array = [NaI_Am_result, NaI_Am15_result, NaI_Am30_result, NaI_Am45_result, NaI_Am60_result, NaI_Am90_result]
    angle_array = [0, 15, 30, 45, 60, 90]
    
    # absolute and intrinsic efficiency against angle
    NaI_abs_eff_angle, NaI_intr_eff_angle = gr_func.angle_efficiencies("NaI", NaI_source_results_array, angle_array)
    
if __name__ == "__main__":
    # there are no arguments to pass
    main()
