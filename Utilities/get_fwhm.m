function fwhm = get_fwhm(t,x)

t_interp = linspace(min(t), max(t), 10000);
y = abs(spline(t, x, t_interp));
peak = 1/2*(min(y) + max(y));
index1 = find(y >= peak, 1, 'first');
index2 = find(y >= peak, 1, 'last');
fwhm = abs(t_interp(index2) - t_interp(index1));


end

