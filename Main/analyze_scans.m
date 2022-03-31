filenames = uigetfile('*.*', 'defname', '.\QW Data\', 'MultiSelect', 'on');

if ~iscell(filenames)
    filenames = {filenames};
end
power_ratio_TE = 13.88;
absorption_per_mm_TE = 0.56;

power_ratio_TM = 11.15;
absorption_per_mm_TM = 0.46;

num_files = length(filenames);
scan_fit = cell(size(filenames));

disp_wavelength_offset_nm = -25;

for ind = 1:num_files
    
    if contains(filenames{ind}, 'TM')
        scan_fit{ind} = scan_analysis_fcn(filenames{ind}, power_ratio_TM, absorption_per_mm_TM, disp_wavelength_offset_nm, 'TM');
        scan_fit{ind}.polarization = 'TM';
    elseif contains(filenames{ind}, 'TE')
        scan_fit{ind} = scan_analysis_fcn(filenames{ind}, power_ratio_TE, absorption_per_mm_TE, disp_wavelength_offset_nm, 'TE');
        scan_fit{ind}.polarization = 'TE';
    end
end
mycellfun = @(scans, field_name) cellfun(@(c) c.(field_name), scans, 'UniformOutput', 0);
TM_indices = find(strcmp(mycellfun(scan_fit, 'polarization'), 'TM'));
TE_indices = find(strcmp(mycellfun(scan_fit, 'polarization'), 'TE'));

scan_fits_TM = scan_fit(TM_indices);
scan_fits_TE = scan_fit(TE_indices);

prompt = 'Save results? [y] or [n]: ';
response = input(prompt, 's');
if strcmp(response, 'y')
    
    if ~isempty(TM_indices)
        save('scan_fits_TM', 'scan_fits_TM');
    end
    if ~isempty(TE_indices)
        save('scan_fits_TE', 'scan_fits_TE');
    end
end