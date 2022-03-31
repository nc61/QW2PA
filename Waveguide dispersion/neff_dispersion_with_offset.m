function [n_eff, beta_1_fs_per_mm, beta_2_fs2_per_mm] = neff_dispersion_with_offset(filename, wavelength_um, wavelength_offset)
c = 3e8;
load(filename)
frequency_input_Hz = c./((wavelength_um - wavelength_offset)*1e-6);
n_eff = interp1(f, neff, frequency_input_Hz, 'spline');
beta_1_fs_per_mm = 1./interp1(f, vg, frequency_input_Hz, 'spline')*1e15/1e3;
beta_2_fs2_per_mm = -(wavelength_um*1e-6).^2./(2*pi*c).*interp1(f, D, frequency_input_Hz, 'spline')*(1e15)^2/1000;

end

