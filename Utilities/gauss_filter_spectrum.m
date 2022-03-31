function smooth_spectrum = gauss_filter_spectrum(wavelength, amplitude)
filter_lambda_nm = -5:wavelength(2) - wavelength(1):5;
gauss_filter = exp(-4*log(2)*filter_lambda_nm.^2./2^2);
smooth_spectrum = conv(amplitude, gauss_filter, 'same');
smooth_spectrum = smooth_spectrum/max(smooth_spectrum);


