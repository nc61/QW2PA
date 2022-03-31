function alpha2_au = TE_2PA_normalized(x1, x2, index_1, index_2, total_thickness_au, kane_energy_au, Eg_GaAs_au, conduction_states, heavy_hole_states, light_hole_states)

kane_energy_gu = kane_energy_au/Eg_GaAs_au;
tangential_energy_gu = @(start_band, start_index, end_band, end_index)(x1 + x2 - (end_band.energies_gu(end_index) - start_band.energies_gu(start_index)));

mu_au = @(start_band, start_index, end_band, end_index) 1./(1./end_band.effective_mass_t_au(end_index) - 1./start_band.effective_mass_t_au(start_index));
kappa_z_0 = @(band, index)sqrt(band.effective_mass_z_au(index)*(band.energies_gu(index) - 1*heaviside(band.energies_gu(index) - 1))*kane_energy_gu);
kappa_t_0 = @(start_band, start_index, end_band, end_index)sqrt(tangential_energy_gu(start_band, start_index, end_band, end_index)*kane_energy_gu.*mu_au(start_band, start_index, end_band, end_index)... 
    .*heaviside(tangential_energy_gu(start_band, start_index, end_band, end_index)));
theta_transition = @(start_band, start_index, end_band, end_index, interband_band, interband_index)acos((1 + kappa_t_0(start_band, start_index, end_band, end_index).^2./kappa_z_0(interband_band, interband_index).^2).^(-1/2)) ... 
    .*heaviside(kappa_t_0(start_band, start_index, end_band, end_index));


theta_lh1 = theta_transition(light_hole_states, 1, conduction_states, 1, conduction_states, 1);
theta_lh2 = theta_transition(light_hole_states, 2, conduction_states, 2, conduction_states, 2);
theta_hh1 = theta_transition(heavy_hole_states, 1, conduction_states, 1, conduction_states, 1);
theta_hh2 = theta_transition(heavy_hole_states, 2, conduction_states, 2, conduction_states, 2);

% theta_lh1 = 0;
% theta_lh2 = 0;
% theta_hh1 = 0;
% theta_hh2 = 0;

lh_m2_average = @(theta)1/96*(17  - 9*cos(2*theta)) + 1/8*(1 - cos(2*theta));
hh_m2_average = @(theta)1/32*(5 + 3*cos(2*theta));


f_scale = 4./(x1.*x2.^2).*(1./x1 + 1./x2).^2;
jdos_energy_c1lh1 = (x1 + x2 - (conduction_states.energies_gu(1) - light_hole_states.energies_gu(1)));
jdos_energy_c2lh2 = (x1 + x2 - (conduction_states.energies_gu(2) - light_hole_states.energies_gu(2)));
jdos_energy_c1hh1 = (x1 + x2 - (conduction_states.energies_gu(1) - heavy_hole_states.energies_gu(1)));
jdos_energy_c2hh2 = (x1 + x2 - (conduction_states.energies_gu(2) - heavy_hole_states.energies_gu(2)));

% figure(1)
% plot(x1 + x2, jdos_energy_c1lh1.*heaviside(jdos_energy_c1lh1), 'r')
% hold on
% plot(x1 + x2, jdos_energy_c2lh2.*heaviside(jdos_energy_c2lh2), 'g')
% plot(x1 + x2, jdos_energy_c1hh1.*heaviside(jdos_energy_c1hh1), 'b')
% plot(x1 + x2, jdos_energy_c2hh2.*heaviside(jdos_energy_c2hh2), 'k')
% hold off

f = f_scale.*(jdos_energy_c1lh1.*heaviside(jdos_energy_c1lh1).*lh_m2_average(theta_lh1) + jdos_energy_c2lh2.*heaviside(jdos_energy_c2lh2).*lh_m2_average(theta_lh2) + ... 
    jdos_energy_c1hh1.*heaviside(jdos_energy_c1hh1).*hh_m2_average(theta_hh1) + jdos_energy_c2hh2.*heaviside(jdos_energy_c2hh2)).*hh_m2_average(theta_hh2);
alpha2_au = f.*K(index_1, index_2, total_thickness_au, kane_energy_au, Eg_GaAs_au);
end


