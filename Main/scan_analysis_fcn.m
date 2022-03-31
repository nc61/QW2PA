
function scan_fit = scan_analysis_fcn(filename, power_ratio, absorption_per_mm, disp_wavelength_offset_nm, polarization)

scan_folder = '..\QW Data';
spectra_folder = '..\QW Data\Spectra';


scan = get_QW_scan(fullfile(scan_folder, filename));
[wavelength_pump_input, amplitude_pump_input] = read_spectrum_csv(fullfile(spectra_folder, scan.pump_spectrum_in_file));
amplitude_pump_input = gauss_filter_spectrum(wavelength_pump_input, amplitude_pump_input);
[wavelength_probe_input, amplitude_probe_input] = read_spectrum_csv(fullfile(spectra_folder, scan.probe_spectrum_in_file));
amplitude_probe_input = gauss_filter_spectrum(wavelength_probe_input, amplitude_probe_input);
[wavelength_pump_output, amplitude_pump_output] = read_spectrum_csv(fullfile(spectra_folder, scan.pump_spectrum_out_file));
amplitude_pump_output = gauss_filter_spectrum(wavelength_pump_output, amplitude_pump_output);
[wavelength_probe_output, amplitude_probe_output] = read_spectrum_csv(fullfile(spectra_folder, scan.probe_spectrum_out_file));
amplitude_probe_output = gauss_filter_spectrum(wavelength_probe_output, amplitude_probe_output);

[wavelength_pump_in_fit, amplitude_pump_in_fit, center_wavelength_pump_in, fwhm_pump_in] = gauss_fit_spectrum(wavelength_pump_input, amplitude_pump_input);
[wavelength_pump_out_fit, amplitude_pump_out_fit, center_wavelength_pump_out, fwhm_pump_out] = gauss_fit_spectrum(wavelength_pump_output, amplitude_pump_output);
[wavelength_probe_in_fit, amplitude_probe_in_fit, center_wavelength_probe_in, fwhm_probe_in] = gauss_fit_spectrum(wavelength_probe_input, amplitude_probe_input);
[wavelength_probe_out_fit, amplitude_probe_out_fit, center_wavelength_probe_out, fwhm_probe_out] = gauss_fit_spectrum(wavelength_probe_output, amplitude_probe_output);

pump_pulsewidth_out_fwhm_fs = get_bwl_pulse(center_wavelength_pump_out, fwhm_pump_out);
probe_pulsewidth_out_fwhm_fs = get_bwl_pulse(center_wavelength_probe_out, fwhm_probe_out);

if strcmp(polarization, 'TE')
[~, pump_beta1_fs_per_mm, pump_dispersion_fs2_per_mm] = neff_dispersion_with_offset('QW_dispersion_3um_ridge_TE_full.mat', 1/1000*center_wavelength_pump_out, disp_wavelength_offset_nm/1000);
[~, probe_beta1_fs_per_mm, probe_dispersion_fs2_per_mm] = neff_dispersion_with_offset('QW_dispersion_3um_ridge_TE_full.mat', 1/1000*center_wavelength_probe_out, disp_wavelength_offset_nm/1000);
elseif strcmp(polarization, 'TM')
[~, pump_beta1_fs_per_mm, pump_dispersion_fs2_per_mm] = neff_dispersion_with_offset('QW_dispersion_3um_ridge_TM_full.mat', 1/1000*center_wavelength_pump_out, disp_wavelength_offset_nm/1000);
[~, probe_beta1_fs_per_mm, probe_dispersion_fs2_per_mm] = neff_dispersion_with_offset('QW_dispersion_3um_ridge_TM_full.mat', 1/1000*center_wavelength_probe_out, disp_wavelength_offset_nm/1000);
end

ssp.delay_fs = scan.delay_fs;
ssp.norm_trans = scan.normalized_transmission;
ssp.pump_pulsewidth_fs_FWHM = 1.1*pump_pulsewidth_out_fwhm_fs;
ssp.probe_pulsewidth_fs_FWHM = probe_pulsewidth_out_fwhm_fs;
ssp.pump_energy_J = 1/82e6*power_ratio*scan.start_power_out_uW*1e-6;
ssp.gvm_fs_per_mm = (pump_beta1_fs_per_mm - probe_beta1_fs_per_mm);
ssp.sample_thickness_mm = 3.6;
ssp.delay_offset_fs = 0;
ssp.gamma_pe_real_minv_Winv = 0;
ssp.gamma_pe_imag_minv_Winv = 0;
ssp.gamma_ee_real_minv_Winv = 0;
ssp.gamma_ee_imag_minv_Winv = 0;
ssp.gamma3_ee_real_minv_W2inv = 0;
ssp.gamma3_ee_imag_minv_W2inv = 0;
ssp.gamma3_pe_real_minv_W2inv = 0;
ssp.gamma3_pe_imag_minv_W2inv = 0.00;

ssp.pump_absorption_per_mm = absorption_per_mm;
ssp.probe_absorption_per_mm = 0;
ssp.pump_gvd_fs2_per_mm = pump_dispersion_fs2_per_mm;
ssp.probe_gvd_fs2_per_mm = probe_dispersion_fs2_per_mm;
ssp.amplitude_offset = 0;
ssp.delay_offset_fs = 0;
ssp.amplitude_offset = 0;

sim_params.num_time_points = 2^9;
sim_params.num_z_steps = 20;

% Probe free parameter
fit_params_guess = [10, 1.2*probe_pulsewidth_out_fwhm_fs, 0];
fit_param_fields = {'gamma_pe_imag_minv_Winv', 'probe_pulsewidth_fs_FWHM', 'delay_offset_fs'};

[~, ~, ~, ~, pump_energy_J, fit_params_out] = fit_scan_split_step(ssp, fit_param_fields, fit_params_guess, sim_params);
gamma_3 = -1*fit_params_out(1);
 
scan_fit.filename = filename;
scan_fit.gamma_probe_free = gamma_3;
scan_fit.pump_energy_J = pump_energy_J;
scan_fit.pump_wavelength_in = center_wavelength_pump_in;
scan_fit.pump_wavelength_out = center_wavelength_pump_out;
scan_fit.probe_wavelength_in = center_wavelength_probe_in;
scan_fit.probe_wavelength_out = center_wavelength_probe_out;
scan_fit.spectral_fwhm_pump_nm = fwhm_pump_out;
scan_fit.spectral_fwhm_pump_nm = fwhm_probe_out;
scan_fit.bwl_pulse_pump_fwhm_fs = pump_pulsewidth_out_fwhm_fs;
scan_fit.probe_pulsewidth_fit_fwhm_fs = fit_params_out(2);
scan_fit.bwl_pulse_probe_fwhm_fs = probe_pulsewidth_out_fwhm_fs;
scan_fit.bwl_pulsewidth_ratio = fit_params_out(2)./probe_pulsewidth_out_fwhm_fs;
scan_fit.sum_wavelength = 1/(1/center_wavelength_probe_out + 1/center_wavelength_pump_out);


 
