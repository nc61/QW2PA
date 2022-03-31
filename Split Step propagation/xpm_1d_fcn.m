function [normalized_probe_energy_out, varargout] = xpm_1d_fcn(t_fs, phi_pump_sqrtW, ...
    phi_probe_sqrtW, sample_thickness_mm, num_z_steps, gamma_ep_minv_Winv, gamma_ee_minv_Winv, gamma3_ep_minv_W2inv, gamma3_ee_minv_W2inv, ...
    gvm_fs_per_mm, ...
    pump_dispersion_fs2_per_mm, probe_dispersion_fs2_per_mm, absorption_pump_per_mm, absorption_probe_per_mm)

if nargout == 1
    output_profiles = 0;
elseif nargout == 4
    output_profiles = 1;
else
    error('Improper number of arguments. Correct numbers should be 1 (normalized energy output only) or 4 (Normalized energy, pump profile, probe_profile, z_mm)')
end

tmax_fs = max(t_fs);
sampling_period_fs = abs(t_fs(2) - t_fs(1));

gamma_ep_mminv_Winv = 1e-3*gamma_ep_minv_Winv;
gamma_ee_mminv_Winv = 1e-3*gamma_ee_minv_Winv;

gamma3_ep_mminv_W2inv = 1e-3*gamma3_ep_minv_W2inv;
gamma3_ee_mminv_W2inv = 1e-3*gamma3_ee_minv_W2inv;

num_time_points = length(t_fs);

% Fourier space for t
omega_prad = 2*pi/sampling_period_fs*[(0:num_time_points/2-1) (-num_time_points/2:-1)]/num_time_points;

% Propagation distance
dz_mm = sample_thickness_mm/num_z_steps;
dt_s = 1e-15*tmax_fs/num_time_points;
initial_probe_energy_W =  dt_s*sum(sum(abs(phi_probe_sqrtW(:,:)).^2));


% Linear propagation operator
D_pump = 1i*(0.5*pump_dispersion_fs2_per_mm*omega_prad.^2 + (gvm_fs_per_mm)*omega_prad + 1i*0.5*absorption_pump_per_mm)*0.5*dz_mm;
D_probe = 1i*(0.5*probe_dispersion_fs2_per_mm*omega_prad.^2 + 0.5*1i*absorption_probe_per_mm)*0.5*dz_mm;

if (output_profiles)
    phi_out_probe_sqrtW = zeros(num_z_steps, num_time_points);
    phi_out_pump_sqrtW = zeros(num_z_steps, num_time_points);
end

z_mm = dz_mm:dz_mm:sample_thickness_mm;


for ind = 1:num_z_steps
    
    lin_step_pump_1_sqrtW = ifft(fft(phi_pump_sqrtW).*exp(D_pump));
    lin_step_probe_1_sqrtW = ifft(fft(phi_probe_sqrtW).*exp(D_probe));
    
    nonlin_step_pump_sqrtW = lin_step_pump_1_sqrtW.*exp((1i*gamma_ee_mminv_Winv*abs(lin_step_pump_1_sqrtW).^2 + 1i*gamma3_ee_mminv_W2inv*abs(lin_step_pump_1_sqrtW).^4)*dz_mm);
    nonlin_step_probe_sqrtW = lin_step_probe_1_sqrtW.*exp((2*1i*gamma_ep_mminv_Winv*abs(lin_step_pump_1_sqrtW).^2 + 3*1i*gamma3_ep_mminv_W2inv*abs(lin_step_pump_1_sqrtW).^4)*dz_mm);
    
    phi_pump_sqrtW = ifft(fft(nonlin_step_pump_sqrtW).*exp(D_pump));
    phi_probe_sqrtW = ifft(fft(nonlin_step_probe_sqrtW).*exp(D_probe));
    
    if (output_profiles)
        phi_out_pump_sqrtW(ind, :) = phi_pump_sqrtW;
        phi_out_probe_sqrtW(ind,:) = phi_probe_sqrtW;
    end
end

final_probe_energy_W =  dt_s*sum(sum(abs(phi_probe_sqrtW(:,:)).^2));
normalized_probe_energy_out = -(initial_probe_energy_W - final_probe_energy_W)/initial_probe_energy_W;

if output_profiles
    varargout{1} = phi_out_pump_sqrtW;
    varargout{2} = phi_out_probe_sqrtW;
    varargout{3} = z_mm;
end

end

