function [x_fit, y_fit, center, fwhm] = gauss_fit_spectrum(x, y)

gauss_fit_fun = @(p, lambda)p(1) + p(2)*exp(-4*log(2)*(lambda - p(3)).^2./p(4).^2);
[ amplitude_guess, center_index_guess] = max(y);
center_wavelength_guess = x(center_index_guess);
offset_guess = min(y);
fwhm_guess = get_fwhm(x, y);

fit_params = lsqcurvefit(gauss_fit_fun, [offset_guess, amplitude_guess, center_wavelength_guess, fwhm_guess], x, y);
[x_fit, y_fit] = get_fit_curve(fit_params, gauss_fit_fun, x, length(x));
center = fit_params(3);
fwhm = fit_params(4);

end

