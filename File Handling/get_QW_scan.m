function scan = get_QW_scan(filename)

text_data = fileread(filename);
lines = regexp(text_data, '\r\n', 'split');

if isempty(lines{end})
    lines = lines(1:end-1);
end

split_lines = regexp(lines, '\t', 'split');

scan.pump_wavelength_nominal_nm = str2double(split_lines{1}{2});
scan.probe_wavelength_nominal_nm = str2double(split_lines{2}{2});
scan.start_power_in_mW = str2double(split_lines{3}{2});
scan.start_power_out_uW = str2double(split_lines{4}{2});
scan.end_power_out_uW = str2double(split_lines{5}{2});
scan.pump_spectrum_in_file = strcat(split_lines{6}{2}, '.csv');
scan.probe_spectrum_in_file = strcat(split_lines{7}{2}, '.csv');
scan.pump_spectrum_out_file = strcat(split_lines{8}{2}, '.csv');
scan.probe_spectrum_out_file = strcat(split_lines{9}{2}, '.csv');

first_col_fun = @(x)x(1);
second_col_fun = @(x)x(2);
data_lines = split_lines(12:end);
first_col = cellfun(first_col_fun, data_lines);
second_col = cellfun(second_col_fun, data_lines);
scan.delay_fs = str2double(first_col);
scan.normalized_transmission = str2double(second_col);

end

