function [alpha2_TM_cm_per_GW, alpha2_TE_cm_per_GW, heavy_hole_states, light_hole_states, conduction_states, Eg_GaAs_au] = AlGaAs_GaAs_QW2PA_fcn(wavelength_pump_um, ... 
    wavelength_probe_um, kane_energy_eV, Ep2_eV, well_width_angstrom, barrier_width_angstrom, barrier_comp, temperature_K, E_increment_meV, z_increment_A, Eg_offset_au, plot_wells)

well_width_au = angstrom_to_au(well_width_angstrom);
barrier_width_au = angstrom_to_au(barrier_width_angstrom);
materials_au = [barrier_width_au 0 barrier_comp; well_width_au 0 0; barrier_width_au 0 barrier_comp];
total_thickness_au = 2*barrier_width_au + well_width_au;
kane_energy_au = eV_to_au(kane_energy_eV);
E_increment_au = meV_to_au(E_increment_meV);
z_increment_au = angstrom_to_au(z_increment_A);
Ep2_au = eV_to_au(Ep2_eV);

if length(wavelength_probe_um) == 1
    wavelength_probe_um = wavelength_probe_um*ones(size(wavelength_pump_um));
end
if length(wavelength_pump_um) == 1
    wavelength_pump_um = wavelength_pump_um*ones(size(wavelength_probe_um));
end

index_1_TM = neff_dispersion('QW_dispersion_3um_ridge_TM.mat', wavelength_probe_um);
index_2_TM = neff_dispersion('QW_dispersion_3um_ridge_TM.mat', wavelength_pump_um);
index_1_TE = neff_dispersion('QW_dispersion_3um_ridge_TE.mat', wavelength_probe_um);
index_2_TE = neff_dispersion('QW_dispersion_3um_ridge_TE.mat', wavelength_pump_um);

[heavy_hole_states, light_hole_states, conduction_states, ~, Eg_GaAs_au] = AlGaAs_states_function(materials_au, temperature_K, E_increment_au, z_increment_au, kane_energy_au, Ep2_au, plot_wells, Eg_offset_au);

x1 = meV_to_au(1240./wavelength_probe_um)/Eg_GaAs_au;
x2 = meV_to_au(1240./wavelength_pump_um)/Eg_GaAs_au;

alpha2_TE_cm_per_GW = au_to_cm_per_GW(TE_2PA_normalized(x1, x2, index_1_TE, index_2_TE, total_thickness_au, kane_energy_au, Eg_GaAs_au, conduction_states, heavy_hole_states, light_hole_states));
alpha2_TM_cm_per_GW = au_to_cm_per_GW(TM_2PA_normalized(x1, x2, index_1_TM, index_2_TM, total_thickness_au, kane_energy_au, Eg_GaAs_au, conduction_states, heavy_hole_states, light_hole_states));

