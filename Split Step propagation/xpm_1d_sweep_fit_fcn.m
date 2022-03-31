function signal = xpm_1d_sweep_fit_fcn(fit_params, fit_param_field_names, split_step_params, simulation_params, delay_fs)
% This function only exists to unwrap the inputs and call the function
% which performs the split step simulations at each delay. Lsqcurvefit
% expects a vector of parameters, specified by the 'initial guess' input.
% This function takes the elements in that vector and maps them to their
% corresponding field in the split_step_params structure. Then the
% split step function is called and the final curve is output to the
% 'signal' variable. Then lsqcurvefit readjusts fit parameters and this
% process repeats until the desired accuracy is achieved.

num_params = length(fit_params);
num_fields = length(fit_param_field_names);

if (num_fields ~= num_params)
    error('You must specify the same number of fields as you do parameters')
end

% Set the split step parameters corresponding to the field names input
% equal to the fit params value. The lsqcurvefit function will vary the
% numbers in these fit params every iteration
for ind = 1:num_params
    split_step_params.(fit_param_field_names{ind}) = fit_params(ind);
end

% This sets the "x" data to be the x-axis data that is passed in
split_step_params.delay_fs = delay_fs;

% Finally, plug these parameters into the simulation which
signal = xpm_1d_sweep_fcn(split_step_params, simulation_params);


end



