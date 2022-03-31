function [alpha2_TM, sum_wavelength_TM, pump_energy_J_TM, alpha2_TE, sum_wavelength_TE, pump_energy_J_TE] = get_QW2PA_data()

%load('TE scans.mat', 'scan_fit');
load('../Final Pump Probe Fits/scan_fits_TE.mat', 'scan_fits_TE');
scans_TE = scan_fits_TE;
num_scans_TE = length(scans_TE);

%load('TM scans.mat', 'scan_fit');
load('../Final Pump Probe Fits/scan_fits_TM.mat', 'scan_fits_TM');
scans_TM = scan_fits_TM;
num_scans_TM = length(scans_TM);

sum_wavelength_TE = zeros(size(scans_TE));
pump_wavelength_TE = zeros(size(scans_TE));
probe_wavelength_TE = zeros(size(scans_TE));
gamma_2pa_only_TE = zeros(size(scans_TE));
gamma_pump_free_TE = zeros(size(scans_TE));
gamma_probe_free_TE = zeros(size(scans_TE));
pump_energy_J_TE = zeros(size(scans_TE));

for ind = 1:num_scans_TE
    sum_wavelength_TE(ind) = scans_TE{ind}.sum_wavelength;
    pump_wavelength_TE(ind) = scans_TE{ind}.pump_wavelength_out;
    probe_wavelength_TE(ind) = scans_TE{ind}.probe_wavelength_out;
    %gamma_2pa_only_TE(ind) = scans_TE{ind}.gamma_2pa_only;
    %gamma_pump_free_TE(ind) = scans_TE{ind}.gamma_pump_free;
    gamma_probe_free_TE(ind) = scans_TE{ind}.gamma_probe_free;
    pump_energy_J_TE(ind) = scans_TE{ind}.pump_energy_J;
end
%gamma_matrix_TE = [gamma_2pa_only_TE; gamma_pump_free_TE; gamma_probe_free_TE];
%gamma_mean_TE = mean(gamma_matrix_TE, 1);

sum_wavelength_TM = zeros(size(scans_TM));
pump_wavelength_TM = zeros(size(scans_TM));
probe_wavelength_TM = zeros(size(scans_TM));
gamma_2pa_only_TM = zeros(size(scans_TM));
gamma_pump_free_TM = zeros(size(scans_TM));
gamma_probe_free_TM = zeros(size(scans_TM));
pump_energy_J_TM = zeros(size(scans_TM));

for ind = 1:num_scans_TM
    sum_wavelength_TM(ind) = scans_TM{ind}.sum_wavelength;
    pump_wavelength_TM(ind) = scans_TM{ind}.pump_wavelength_out;
    probe_wavelength_TM(ind) = scans_TM{ind}.probe_wavelength_out;
    %gamma_2pa_only_TM(ind) = scans_TM{ind}.gamma_2pa_only;
    %gamma_pump_free_TM(ind) = scans_TM{ind}.gamma_pump_free;
    gamma_probe_free_TM(ind) = scans_TM{ind}.gamma_probe_free;
    pump_energy_J_TM(ind) = scans_TM{ind}.pump_energy_J;
end
%gamma_matrix_TM = [gamma_2pa_only_TM; gamma_pump_free_TM; gamma_probe_free_TM];
%gamma_mean_TM = mean(gamma_matrix_TM, 1);

alpha2_TE = 1/10*gamma_probe_free_TE*2.*sqrt(get_mode_area(probe_wavelength_TE, 'TE')).*sqrt(get_mode_area(pump_wavelength_TE, 'TE'));
alpha2_TM = 1/10*gamma_probe_free_TM*2.*sqrt(get_mode_area(probe_wavelength_TM, 'TM')).*sqrt(get_mode_area(pump_wavelength_TM, 'TM'));