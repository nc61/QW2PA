function alpha2_au = TM_2PA_normalized(x1, x2, index_1, index_2, total_thickness_au, kane_energy_au, Eg_GaAs_au, conduction_states, heavy_hole_states, light_hole_states)

kane_energy_gu = kane_energy_au/Eg_GaAs_au;
kane_parameter_au = sqrt(kane_energy_au/2);

for ind = 1:size(conduction_states.wavefunctions, 1)
    conduction_states.wavefunction_derivative(ind,:) = [0 diff(conduction_states.wavefunctions(ind,:))]; 
end
for ind = 1:size(heavy_hole_states.wavefunctions, 1)
    heavy_hole_states.wavefunction_derivative(ind,:) = [0 diff(heavy_hole_states.wavefunctions(ind,:))]; 
end
for ind = 1:size(light_hole_states.wavefunctions, 1)
    light_hole_states.wavefunction_derivative(ind,:) = [0 diff(light_hole_states.wavefunctions(ind,:))]; 
end

intersubband_element = @(band, subband_start, subband_end)1/kane_parameter_au*(1/2*(trapz((1./band.effective_mass_z_au(subband_start,:) - band.F).*band.wavefunction_derivative(subband_start,:).*band.wavefunctions(subband_end, :)) ...
    - trapz((1./band.effective_mass_z_au(subband_end,:) - band.F).*band.wavefunction_derivative(subband_end,:).*band.wavefunctions(subband_start,:))));

% intersubband_element = @(band, subband_start, subband_end)1/kane_parameter_au*(1/2*(trapz((1./band.effective_mass_z_au(subband_start,:)).*band.wavefunction_derivative(subband_start,:).*band.wavefunctions(subband_end, :)) ...
%     - trapz((1./band.effective_mass_z_au(subband_end,:)).*band.wavefunction_derivative(subband_end,:).*band.wavefunctions(subband_start,:))));

energy_dispersion = @(band, index, kappa)band.energies_gu(index) + kappa.^2./(band.effective_mass_t_au(index).*kane_energy_gu);

tangential_energy_gu = @(start_band, start_index, end_band, end_index)(x1 + x2 - (end_band.energies_gu(end_index) - start_band.energies_gu(start_index)));

mu_t_au = @(start_band, start_index, end_band, end_index) 1./(1./end_band.effective_mass_t_au(end_index) - 1./start_band.effective_mass_t_au(start_index));

kappa_t_0 = @(start_band, start_index, end_band, end_index)sqrt(tangential_energy_gu(start_band, start_index, end_band, end_index)*kane_energy_gu.*mu_t_au(start_band, start_index, end_band, end_index)... 
    .*heaviside(tangential_energy_gu(start_band, start_index, end_band, end_index)));

kappa_z_0 = @(band, index)sqrt(band.effective_mass_z_au(index)*(band.energies_gu(index) - 1*heaviside(band.energies_gu(index) - 1))*kane_energy_gu);

theta_transition = @(start_band, start_index, end_band, end_index, interband_band, interband_index)acos((1 + kappa_t_0(start_band, start_index, end_band, end_index).^2./kappa_z_0(interband_band, interband_index).^2).^(-1/2)) ... 
    .*heaviside(kappa_t_0(start_band, start_index, end_band, end_index));

detuning = @(start_band, start_index, end_band, end_index, int_band, int_index, x)energy_dispersion(int_band, int_index, kappa_t_0(start_band, start_index, end_band, end_index)) -  ... 
    energy_dispersion(start_band, start_index, kappa_t_0(start_band, start_index, end_band, end_index)) - x;
detuning_denominator = @(start_band, start_index, end_band, end_index, int_band, int_index) ...
    1./detuning(start_band, start_index, end_band, end_index, int_band, int_index, x1) + 1./detuning(start_band, start_index, end_band, end_index, int_band, int_index, x2);


S_LH1C2_a = S_function(light_hole_states, 1, 2, (2/3), @cos);


S_LH1C2_b = S_function(light_hole_states, 1, 2, (1/6), @sin);


S_LH2C1_a = S_function(light_hole_states, 2, 1, (2/3), @cos);
S_LH2C1_b = S_function(light_hole_states, 2, 1, (1/6), @sin);

S_HH1C2 = S_function(heavy_hole_states, 1, 2, (1/2), @sin);
S_HH2C1 = S_function(heavy_hole_states, 2, 1, (1/2), @sin);

S_tot = S_LH1C2_a + S_LH1C2_b + S_LH2C1_a + S_LH2C1_b + S_HH1C2 + S_HH2C1;

f = 1./(x1.*x2.^2).*S_tot*kane_energy_gu;
alpha2_au = f.*K(index_1, index_2, total_thickness_au, kane_energy_au, Eg_GaAs_au);


function [S, transition_1, transition_2] = S_function(start_band, start_index, end_index, interband_factor, interband_theta_function)
% Theta value for when interband transition occurs from H1 -> C1
theta_1 = theta_transition(start_band, start_index, conduction_states, end_index, conduction_states, start_index);

%Theta value for when interband transition occurs from H2 -> C2
theta_2 = theta_transition(start_band, start_index, conduction_states, end_index, conduction_states, end_index);

% Interband transition to conduction band followed by intersubband
% transition in conduction band
transition_1 = interband_theta_function(theta_1).*intersubband_element(conduction_states, start_index, end_index).*detuning_denominator(start_band, start_index, conduction_states, end_index, conduction_states, start_index) ...
    .*heaviside(tangential_energy_gu(start_band,start_index,conduction_states,end_index));

%Intersubband transition within valence band followed by interband
%transition
transition_2 = interband_theta_function(theta_2).*intersubband_element(start_band, start_index, end_index).*detuning_denominator(start_band, start_index, conduction_states, end_index, start_band, end_index) ...
    .*heaviside(tangential_energy_gu(start_band,start_index,conduction_states,end_index));



S = mu_t_au(start_band, start_index, conduction_states, end_index).*interband_factor.*abs(transition_1 + transition_2).^2;
%S = interband_factor.*abs(transition_1 + transition_2).^2;
% 
% figure(2)
% plot(x1 + x2, transition_1, 'b')
% hold on
% plot(x1 + x2, transition_2, 'r')
% plot(x1 + x2, S, 'k')
% hold off

end

end




