function pulsewidth_fwhm_fs = get_bwl_pulse(center_nm, fwhm_nm)

c = 3e8;
pulsewidth_fwhm_fs = 1e6*0.44*center_nm^2/(c*fwhm_nm);

end

