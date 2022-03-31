function [signal, varargout] = xpm_1d_sweep_fcn(ssp, sim_params)
% ssp = split step parameters. Name shortened for easier accessing of
% fields throughout the code.
% ssp must have the following fields:
%   pump_pulsewidth_fs_FWHM
%   probe_pulsewidth_fs_FWHM
%   pump_energy_J
%   gvm_fs_per_mm
%   sample_thickness_mm
%   delay_offset_fs
%   gamma_pe_mminv_Winv
%   gamma_ee_mminv_Winv
%   pump_absorption_per_mm
%   probe_absorption_per_mm
%   pump_gvd_fs2_per_mm
%   probe_gvd_fs2_per_mm
%
%   sim_params is the structure containing information about how the
%   simulation is performed. The fields are:
%
%   num_time_points: Should be a power of 2. This is the number of points
%       used to sample the signal in the time domain
%   num_z_steps: How many steps in the z direction to split the propagation
%       into

if nargout == 6
    store_profiles = 1;
else
    store_profiles = 0;
end

probe_pulsewidth_fs = ssp.probe_pulsewidth_fs_FWHM/(2*sqrt(log(2)));
pump_pulsewidth_fs = ssp.pump_pulsewidth_fs_FWHM/(2*sqrt(log(2)));

E_probe = 1;

pump_peak_power_W = ssp.pump_energy_J/(pump_pulsewidth_fs*1e-15*sqrt(pi));
probe_peak_power_W = E_probe/(probe_pulsewidth_fs*1e-15*sqrt(pi));

time_shift_fs = -ssp.gvm_fs_per_mm*ssp.sample_thickness_mm;

delay_fs = ssp.delay_fs - ssp.delay_offset_fs;
signal = zeros(size(delay_fs));

gamma_pe_minv_Winv = (ssp.gamma_pe_real_minv_Winv - 1i*ssp.gamma_pe_imag_minv_Winv);
gamma_ee_minv_Winv = (ssp.gamma_ee_real_minv_Winv - 1i*ssp.gamma_ee_imag_minv_Winv);
gamma3_pe_minv_W2inv = (ssp.gamma3_pe_real_minv_W2inv - 1i*ssp.gamma3_pe_imag_minv_W2inv);
gamma3_ee_minv_W2inv = (ssp.gamma3_ee_real_minv_W2inv - 1i*ssp.gamma3_ee_imag_minv_W2inv);

if (store_profiles)
    phi_out_pump = cell(1,length(delay_fs));
    phi_out_probe = cell(1,length(delay_fs));
    z_mm = cell(1, length(delay_fs));
    t_fs_cellarray = cell(1, length(delay_fs));
end

for ind = 1:length(delay_fs)
    
    % Number of fourier modes
    tmax_fs = 20*max([probe_pulsewidth_fs, pump_pulsewidth_fs]) + (2*abs(time_shift_fs) + abs(delay_fs(ind)));
    
    % points for x and y
    t_fs = (tmax_fs/sim_params.num_time_points)*(-sim_params.num_time_points/2+1:sim_params.num_time_points/2);
    
    phi_pump = sqrt(pump_peak_power_W).*exp(-(t_fs - delay_fs(ind)).^2./(2*pump_pulsewidth_fs.^2));
    phi_probe = sqrt(probe_peak_power_W).*exp(-(t_fs).^2./(2*probe_pulsewidth_fs.^2));
    
    if (store_profiles)
        [normalized_probe_energy_out, phi_out_pump{ind}, phi_out_probe{ind}, z_mm{ind}] = xpm_1d_fcn(t_fs, phi_pump, phi_probe, ssp.sample_thickness_mm, sim_params.num_z_steps, gamma_pe_minv_Winv, gamma3_ee_minv_W2inv, gamma3_pe_minv_W2inv, ...
            gamma_ee_minv_Winv, ssp.gvm_fs_per_mm, ssp.pump_gvd_fs2_per_mm, ssp.probe_gvd_fs2_per_mm, ssp.pump_absorption_per_mm, ssp.probe_absorption_per_mm);
    else
        [normalized_probe_energy_out] = xpm_1d_fcn(t_fs, phi_pump, phi_probe, ssp.sample_thickness_mm, sim_params.num_z_steps, gamma_pe_minv_Winv, gamma3_ee_minv_W2inv, gamma3_pe_minv_W2inv, ...
            gamma_ee_minv_Winv, ssp.gvm_fs_per_mm, ssp.pump_gvd_fs2_per_mm, ssp.probe_gvd_fs2_per_mm, ssp.pump_absorption_per_mm, ssp.probe_absorption_per_mm);
    end
    
    signal(ind) = normalized_probe_energy_out;
    
end

signal = 1 + signal + ssp.amplitude_offset;

if (store_profiles)
    varargout{1} = phi_out_pump;
    varargout{2} = phi_out_probe;
    varargout{3} = z_mm;
    varargout{4} = delay_fs;
    varargout{5} = t_fs_cellarray;
end
end

