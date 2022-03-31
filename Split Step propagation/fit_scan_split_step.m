function [scan_delay_fs, scan_norm_trans, delay_fit_fs, norm_trans_fit, pump_energy_J, fit_params_out] = fit_scan_split_step(split_step_params, fit_param_fields, fit_params_guess, sim_params)

% Since the offset parameters are not known a priori, they are only fitting
% parameters. Therefore they haven't been defined as fields of the
% structure so we want to make sure they are initialized. Otherwise the
% next line where we check field names will fail.
split_step_params.delay_offset_fs = 0;
split_step_params.amplitude_offset = 0;

% Verify that the fields chosen as fitting parameters are actually fields
% of the split step parameters struct
check_field_names(split_step_params, fit_param_fields)

% Define the anonymous function to use for least squares fitting
fit_fun = @(fit_params, x) xpm_1d_sweep_fit_fcn(fit_params, fit_param_fields, split_step_params, sim_params, x);

% Get the fitting parameters whose simulation results most closely match
% the passed in data (delay_fs and norm_trans)
fit_params_out = lsqcurvefit(fit_fun, fit_params_guess, split_step_params.delay_fs, split_step_params.norm_trans);

% Use the fitting results to generate a curve which can be plotted on top
% of the original data
[delay_fit_fs, norm_trans_fit] = get_fit_curve(fit_params_out, fit_fun, split_step_params.delay_fs, 500);

% Readjust the scan delay output so it is centered at zero delay
scan_delay_fs = split_step_params.delay_fs;

% Readjust the normalized transmission so the baseline is at 100%
scan_norm_trans = split_step_params.norm_trans;

pump_energy_J = split_step_params.pump_energy_J;

end

