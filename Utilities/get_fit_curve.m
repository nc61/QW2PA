function [x, y] = get_fit_curve(fit_params, fitting_function, input_xdata, num_points)

x = linspace(min(input_xdata), max(input_xdata), num_points);
y = fitting_function(fit_params, x);


end

