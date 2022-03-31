function [wavelength, amplitude] = read_spectrum_csv(filename)

headerlines = csvread(filename, 2, 0, [2,0,2,0]);
data = csvread(filename, headerlines + 4);
wavelength = data(:,1);
amplitude = data(:,2);