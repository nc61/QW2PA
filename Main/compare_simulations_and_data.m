sum_wavelength_um = linspace(0.5,0.81, 1000);
wavelength_pump_um = 1.96*ones(size(sum_wavelength_um));
wavelength_pump_um_deg = sum_wavelength_um*2;

kane_energy_14_band_eV = 28.9;
Ep2_14_band_eV = 6;

kane_energy_14_band_kane_free_eV = 29.8;
Ep2_14_band_kane_free_eV = 8;

kane_energy_14_band_wrong_eV = 25.7;
Ep2_14_band_wrong_eV = 6;

kane_energy_8_band_eV = 28.9;
Ep2_8_band_eV = 0;

well_width_angstrom_nominal = 80;
barrier_width_angstrom_nominal = 60;

ratio_g = 0.98;
ratio_a = 1.03;

growth_rate_Ga_nominal = 0.745*2.83;
growth_rate_Al_nominal = 0.350*2.83;
growth_rate_Ga = ratio_g*growth_rate_Ga_nominal;
growth_rate_Al = ratio_a*growth_rate_Al_nominal;
barrier_comp = growth_rate_Al/(growth_rate_Al+growth_rate_Ga);
barrier_time = barrier_width_angstrom_nominal/(growth_rate_Al_nominal + growth_rate_Ga_nominal);
well_time = well_width_angstrom_nominal/growth_rate_Ga_nominal;
well_width_angstrom = well_time*growth_rate_Ga;
barrier_width_angstrom = barrier_time*(growth_rate_Ga + growth_rate_Al);
total_width = barrier_width_angstrom + well_width_angstrom;

plot_wells = 1;

temperature_K = 295;
E_increment_meV = 0.1;
z_increment_A = 0.1;
Eg_offset_meV = 0;
wavelength_shift_nm = 0;
Eg_offset_au = meV_to_au(Eg_offset_meV);

final_scaling = 1.0;

wavelength_probe_um = 1./(1./sum_wavelength_um - 1./wavelength_pump_um);

[alpha2_TM_cm_per_GW_sim_14_band, alpha2_TE_cm_per_GW_sim_14_band, heavy_hole_states_14_band, light_hole_states_14_band, conduction_states_14_band, Eg_GaAs_au] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, wavelength_probe_um, kane_energy_14_band_eV, Ep2_14_band_eV, well_width_angstrom, barrier_width_angstrom, ... 
    barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells);

[alpha2_TM_cm_per_GW_sim_14_band_no_correction, ~, ~, ~, ~, ~] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, wavelength_probe_um, kane_energy_14_band_eV, Ep2_14_band_eV, well_width_angstrom_nominal, barrier_width_angstrom_nominal, ... 
    barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells);

[alpha2_TM_cm_per_GW_sim_14_band_kane_free, alpha2_TE_cm_per_GW_sim_14_band_kane_free, heavy_hole_states, light_hole_states, conduction_states, Eg_GaAs_au] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, wavelength_probe_um, kane_energy_14_band_kane_free_eV, Ep2_14_band_kane_free_eV, well_width_angstrom, barrier_width_angstrom, ... 
    barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells);

[alpha2_TM_cm_per_GW_sim_14_band_wrong, alpha2_TE_cm_per_GW_sim_14_band_wrong, heavy_hole_states_wrong, light_hole_states_wrong, conduction_states_wrong, Eg_GaAs_au] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, wavelength_probe_um, kane_energy_14_band_wrong_eV, Ep2_14_band_wrong_eV, well_width_angstrom, barrier_width_angstrom, ... 
    barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells);

[alpha2_TM_cm_per_GW_sim_8_band, alpha2_TE_cm_per_GW_sim_8_band, heavy_hole_states_8_band, light_hole_states_8_band, conduction_states_8_band, Eg_GaAs_au] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, wavelength_probe_um, kane_energy_8_band_eV, Ep2_8_band_eV, well_width_angstrom, barrier_width_angstrom, ... 
    barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells);

[alpha2_TM_cm_per_GW_data, sum_wavelength_TM_nm_data, pump_energy_J_TM, alpha2_TE_cm_per_GW_data, sum_wavelength_TE_nm_data, pump_energy_J_TE] = get_QW2PA_data();

sum_wavelength_nm_sim_shifted = 1000*sum_wavelength_um - wavelength_shift_nm;

dlambda = sum_wavelength_nm_sim_shifted(2) - sum_wavelength_nm_sim_shifted(1);
filter_width = 4;
filter_wavelength = -3*filter_width:dlambda:3*filter_width;
gaussian_blur = 1/sqrt(pi)/filter_width*exp(-(filter_wavelength).^2/(filter_width)^2);

alpha2_TM_sim_14_band = alpha2_TM_cm_per_GW_sim_14_band;
alpha2_TM_sim_14_band_no_correction = alpha2_TM_cm_per_GW_sim_14_band_no_correction;
alpha2_TM_sim_14_band_kane_free = alpha2_TM_cm_per_GW_sim_14_band_kane_free;
alpha2_TM_sim_14_band_wrong = alpha2_TM_cm_per_GW_sim_14_band_wrong;
alpha2_TM_sim_8_band = alpha2_TM_cm_per_GW_sim_8_band;

alpha2_TM_sim_14_band_blurred = dlambda*conv(alpha2_TM_sim_14_band, gaussian_blur, 'same');
alpha2_TM_sim_14_band_no_correction_blurred = dlambda*conv(alpha2_TM_sim_14_band_no_correction, gaussian_blur, 'same');
alpha2_TM_sim_14_band_kane_free_blurred = dlambda*conv(alpha2_TM_sim_14_band_kane_free, gaussian_blur, 'same');
alpha2_TM_sim_8_band_blurred = dlambda*conv(alpha2_TM_sim_8_band, gaussian_blur, 'same');

[sum_wavelength_TM_high_energy, alpha2_TM_high_energy, sum_wavelength_TM_low_energy, alpha2_TM_low_energy, ... 
    sum_wavelength_TE_high_energy, alpha2_TE_high_energy, sum_wavelength_TE_low_energy, alpha2_TE_low_energy] = ... 
    split_energy_alpha2(pump_energy_J_TM, pump_energy_J_TE, alpha2_TM_cm_per_GW_data, alpha2_TE_cm_per_GW_data, ... 
    sum_wavelength_TM_nm_data, sum_wavelength_TE_nm_data);

x0 = 0;
y0 = 0;
width = 3.375;
height = 2.53;

figure('Units', 'inches', 'Position', [x0 y0 width height], 'PaperPositionMode', 'auto')
plot(sum_wavelength_nm_sim_shifted, alpha2_TM_sim_14_band, 'k', 'linewidth', 1.3)
hold on
plot(sum_wavelength_nm_sim_shifted, alpha2_TM_sim_14_band_blurred, 'b--', 'linewidth', 1.0)
plot(sum_wavelength_nm_sim_shifted, alpha2_TM_sim_14_band_no_correction, 'm:', 'linewidth', 1.8)
plot(sum_wavelength_nm_sim_shifted, alpha2_TM_sim_14_band_no_correction_blurred, 'm--', 'linewidth', 1.3)

plot(sum_wavelength_TM_high_energy, alpha2_TM_high_energy, 'k^')
plot(sum_wavelength_TM_low_energy, alpha2_TM_low_energy, 'ko')
hold off
xlim([732 796])

ylabel('$\alpha_2$ (cm/GW)', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');
xlabel('Sum Wavelength (nm)', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');
title('TM-TM', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');

figure('Units', 'inches', 'Position', [x0 y0 width height], 'PaperPositionMode', 'auto')
plot(sum_wavelength_TE_high_energy, alpha2_TE_high_energy, 'k^')
hold on
plot(sum_wavelength_TE_low_energy, alpha2_TE_low_energy, 'ko')
plot(sum_wavelength_nm_sim_shifted, 5*alpha2_TE_cm_per_GW_sim_14_band, 'k', 'linewidth', 2)
hold off
xlim([732 796])

ylabel('$\alpha_2$ (cm/GW)', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');
xlabel('Sum Wavelength (nm)', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');
title('TE-TE', 'interpreter', 'latex', 'fontUnits', 'points', 'FontSize', 10, 'FontName', 'Times');

heavy_hole_z_inv_mass = -(6.98 - 2*2.08);

F = conduction_states.F + 1;
gamma1 = 1/2*(-heavy_hole_z_inv_mass + 1 + light_hole_states.F);
gamma2 = 1/2*(gamma1 - (-heavy_hole_z_inv_mass));

function [sum_wavelength_TM_high_energy, alpha2_TM_high_energy, sum_wavelength_TM_low_energy, alpha2_TM_low_energy, ... 
    sum_wavelength_TE_high_energy, alpha2_TE_high_energy, sum_wavelength_TE_low_energy, alpha2_TE_low_energy] = ... 
    split_energy_alpha2(pump_energy_J_TM, pump_energy_J_TE, alpha2_TM, alpha2_TE, sum_wavelength_TM, sum_wavelength_TE)
    
    low_energy_condition_TM = abs(pump_energy_J_TM*1e12 - 3.5) < 0.5;
    alpha2_TM_low_energy = alpha2_TM(low_energy_condition_TM);
    sum_wavelength_TM_low_energy = sum_wavelength_TM(low_energy_condition_TM);
    
    low_energy_condition_TE = abs(pump_energy_J_TE*1e12 - 4.2) < 0.5;
    alpha2_TE_low_energy = alpha2_TE(low_energy_condition_TE);
    sum_wavelength_TE_low_energy = sum_wavelength_TE(low_energy_condition_TE);
    
    high_energy_condition_TM = abs(pump_energy_J_TM*1e12 - 4.9) < 0.5;
    alpha2_TM_high_energy = alpha2_TM(high_energy_condition_TM);
    sum_wavelength_TM_high_energy = sum_wavelength_TM(high_energy_condition_TM);
    
    high_energy_condition_TE = abs(pump_energy_J_TE*1e12 - 6) < 0.5;
    alpha2_TE_high_energy = alpha2_TE(high_energy_condition_TE);
    sum_wavelength_TE_high_energy = sum_wavelength_TE(high_energy_condition_TE);
end
