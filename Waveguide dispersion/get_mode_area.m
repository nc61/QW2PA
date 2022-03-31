function area_out = get_mode_area(wavelength_nm, polarization)

areas_TE_um2 = [3.0728 3.08758 3.12199 3.15636 3.1711];
areas_TM_um2 = [3.12146 3.13809 3.17695 3.21596 3.23275];
wavelength = [1150 1178 1250 1320 1350];

area_out = -1;

if strcmpi(polarization, 'TE')
    
    if (wavelength_nm < 1800)
        area_out = interp1(wavelength, areas_TE_um2, wavelength_nm, 'spline');
    else
        area_out = 3.47631;
    end
    
elseif strcmpi(polarization, 'TM')
    
    if (wavelength_nm < 1800)
        area_out = interp1(wavelength, areas_TM_um2, wavelength_nm, 'spline');
        
    else
        area_out = 3.60;
    end
end


