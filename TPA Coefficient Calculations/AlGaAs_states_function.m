function [heavy_hole_states, light_hole_states, conduction_states, z_au, Eg_GaAs_au] = AlGaAs_states_function(materials_au, temperature_K, E_increment_au, z_increment_au, kane_energy_au, Ep2_au, plot_wells, Eg_offset_au)
% For a given matrix of materials, calculate the bound states in the
% conduction and valence bands.
% 
%   calls: schrodinger.m, solver for schrodinger equation
%   
%   materials: each row of the materials matrix is a different material
%   with thickness in Angstroms. The first column is the composition (Ex:
%   0.32) and the third column is the thickness. The second column is
%   reserved for the future when this supports more material systems
%   plot_wells: Display the wavefunctions for visualization
%   Eincrement: Passed through to schrodinger.m, energy resolution of
%   schrodinger solver (very important)


% Go through the materials matrix and find the composition at every
% individual Angstrom step in z

total_length_au = sum(materials_au(:,1));
z_au = 1:z_increment_au:total_length_au;

current_position_au = 0;
x = -1*ones(size(z_au));
for ind = 1:size(materials_au, 1)
    layer_thickness_au = materials_au(ind, 1);
    x(z_au >= current_position_au & z_au <= current_position_au + layer_thickness_au) = materials_au(ind, 3);
    current_position_au = current_position_au + layer_thickness_au;
end

if sum(x == -1) == 1
    error("Layers not fully defined")
end


% Bandgaps (au)
alpha = 8.95e-1;
beta = 538;
Eg_GaAs_au = meV_to_au(1000*1.519 - alpha*temperature_K^2/(temperature_K + beta)) + Eg_offset_au;

alpha = 4.0e-1;
beta = 241;
Eg2_GaAs_au = meV_to_au(1000*4.509 - alpha*temperature_K^2/(temperature_K + beta));


% Effective masses
% Adachi properties of Aluminum Gallium Arsenide
% m_e_eff_AlGaAs_z = (0.0665 + 0.083.*x).*m0;

gamma_1 = @(x) 6.98 + (3.76 - 6.98).*x;
gamma_2 = @(x) 2.06 + (0.82 - 2.06).*x;

% V (conduction band offset)
Q_v = 0.33;

Ec_offset_au = meV_to_au(1395.*x.*(1 - Q_v));
Ev_offset_au = meV_to_au(-1395.*x.*Q_v);

% Applied voltage, should be added as a parameter and not automatically set
% to 0
V_app_au = 0;

% Potential barriers with adjustments based on applied bias.
V_c_au = Ec_offset_au - V_app_au*z_au/max(z_au);
V_hh_au = Ev_offset_au - V_app_au*z_au/max(z_au);
V_lh_au = Ev_offset_au - V_app_au*z_au/max(z_au);
V_c2_au = Ec_offset_au*20;

m_inv_c_literature = 1/0.0635;
m_inv_lh_literature = -1/0.077;

Delta_au = meV_to_au(0.341*1000);
Delta_2_au = eV_to_au(0.171);


m_inv_c_eff_0 = 1 + kane_energy_au/3.*(2./Eg_GaAs_au + 1./(Eg_GaAs_au + Delta_au));
m_inv_c_eff_0 = 1 + kane_energy_au/3.*(2./Eg_GaAs_au + 1./(Eg_GaAs_au + Delta_au)) - Ep2_au/3.*(2./(-Eg_GaAs_au + Delta_2_au + Eg2_GaAs_au) + 1./(-Eg_GaAs_au + Eg2_GaAs_au));

m_inv_lh_eff_0 = 1 + 2*kane_energy_au/3.*(1./(-Eg_GaAs_au));

F_c = (- m_inv_c_eff_0 + m_inv_c_literature);
F_lh = (- m_inv_lh_eff_0 + m_inv_lh_literature);


m_c_eff_AlGaAs_z_au_fcn = @(E_au)1./(1 + F_c + kane_energy_au/3.*(2./(E_au + Eg_GaAs_au - V_lh_au) + 1./(E_au + Eg_GaAs_au + Delta_au - V_lh_au)));
m_c_eff_AlGaAs_z_au_fcn = @(E_au)1./(1 + F_c + kane_energy_au/3.*(2./(E_au + Eg_GaAs_au - V_lh_au) + 1./(E_au + Eg_GaAs_au + Delta_au - V_lh_au))  -  ...
    Ep2_au/3.*(2./(-E_au + Delta_2_au + Eg2_GaAs_au + V_c_au) + 1./(-E_au + Eg2_GaAs_au + V_c_au)));

m_lh_eff_AlGaAs_z_au_fcn = @(E_au)1./(1 + F_lh + 2*kane_energy_au/3.*(1./(E_au - Eg_GaAs_au - V_c_au)));
m_hh_eff_AlGaAs_z_au_fcn = @(E_au)ones(size(E_au)).*-1./(gamma_1(x) - 2*gamma_2(x));

% m_c_eff_AlGaAs_z_au_fcn = @(E_au)(E_au./E_au).*(0.063 + 0.083*x);
% m_lh_eff_AlGaAs_z_au_fcn = @(E_au)(E_au./E_au).*-1./(gamma_1(x) + 2*gamma_2(x));
% m_hh_eff_AlGaAs_z_au_fcn = @(E_au)(E_au./E_au).*-1./(gamma_1(x) - 2*gamma_2(x));



% Use the schrodinger solver to solve for electron, light hole and heavy
% hole states
[wavefunctions_lh, energies_lh_au, dz_lh_au] = schrodinger_solver_variable_meff_au(V_lh_au, z_au, m_lh_eff_AlGaAs_z_au_fcn, E_increment_au);
[wavefunctions_hh, energies_hh_au, dz_hh_au] = schrodinger_solver_variable_meff_au(V_lh_au, z_au, m_hh_eff_AlGaAs_z_au_fcn, E_increment_au);
[wavefunctions_c, energies_c_au, dz_c_au] = schrodinger_solver_variable_meff_au(V_c_au, z_au, m_c_eff_AlGaAs_z_au_fcn, E_increment_au);

% Shift the conduction band energies upward by the bandgap. Now the top of
% the valence bands are at 0 energy and the bottom of the conduction band
% is at 1424meV (GaAs bandgap). THIS SHOULD BE CHANGED IF THE WELLS AREN'T
% MADE OUT OF GAAS
energies_c_au = energies_c_au + Eg_GaAs_au;
V_c_au = V_c_au + Eg_GaAs_au;

% Save the heavy hole states into a structure for output
heavy_hole_states.wavefunctions = wavefunctions_hh;
heavy_hole_states.energies_gu = energies_hh_au/Eg_GaAs_au;
heavy_hole_states.dz_au = dz_hh_au;
heavy_hole_states.effective_mass_z_au = m_hh_eff_AlGaAs_z_au_fcn(energies_hh_au);
heavy_hole_states.effective_mass_t_au = (energies_hh_au./energies_hh_au).*-1./(gamma_1(0) + gamma_2(0));
heavy_hole_states.F = -1 + m_hh_eff_AlGaAs_z_au_fcn(0);

% Save the light hole states into a structure for output
light_hole_states.wavefunctions = wavefunctions_lh;
light_hole_states.energies_gu = energies_lh_au/Eg_GaAs_au;
light_hole_states.dz_au = dz_lh_au;
light_hole_states.effective_mass_z_au = m_lh_eff_AlGaAs_z_au_fcn(energies_lh_au);
light_hole_states.effective_mass_t_au = (energies_lh_au./energies_lh_au).*-1./(gamma_1(0) - gamma_2(0));
light_hole_states.F = F_lh;

% Save the electron states into a structure for output
conduction_states.wavefunctions = wavefunctions_c;
conduction_states.energies_gu = energies_c_au/Eg_GaAs_au;
conduction_states.dz_au = dz_c_au;
conduction_states.effective_mass_z_au = m_c_eff_AlGaAs_z_au_fcn(energies_c_au - Eg_GaAs_au);
conduction_states.effective_mass_t_au = (energies_c_au./energies_c_au).*0.0635;
conduction_states.F = F_c;

% Plot the output if desired
if ~isempty(energies_hh_au) && (plot_wells ~= 0) 
    
    % Kind of a workaround to make sure we don't cover up this plot with
    % whavever we are plotting afterward
    figure(plot_wells)

    scale = .5*energies_hh_au(1)/max(wavefunctions_hh(1,:));
    for ind = 1:size(wavefunctions_hh, 1)
        %diagram = ones(size(z)).*energies_hh(ind);
        diagram = energies_hh_au(ind) + scale*wavefunctions_hh(ind, :);
        plot(z_au*1e9, diagram, 'b');
        hold on;
    end
    
    scale = .5*energies_lh_au(1)/max(wavefunctions_lh(1,:));
    for ind = 1:size(wavefunctions_lh, 1)
        %diagram = ones(size(z)).*energies_lh(ind);
        diagram = energies_lh_au(ind) + scale*wavefunctions_lh(ind, :);
        plot(z_au*1e9, diagram, 'r');
        hold on;
    end
    
    scale = .5*((energies_c_au(1) - Eg_GaAs_au)/max(wavefunctions_c(1,:)));
    for ind = 1:size(wavefunctions_c, 1)
        %diagram = ones(size(z)).*energies_c(ind);
        diagram = energies_c_au(ind) + scale*wavefunctions_c(ind, :);
        plot(z_au*1e9, diagram, 'g');
        hold on;
    end
    
    plot(z_au*1e9, V_c_au, 'k--');
    title('Solutions to Schrodinger equation for GaAs/AlGaAs structures'), xlabel('z (nm)'), ylabel('E (au)');
    plot(z_au*1e9, V_hh_au, 'k--');
    hold off;
end


end

