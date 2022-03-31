function scan = get_pump_probe_scan(directory_sheet, scan_name, micrometer_units)

scan = read_scan_directory(directory_sheet, scan_name);

[scan.delay_fs, scan.norm_trans] = read_scan('filename', scan_name, 'plot', 'off', 'normalizeToProbe', 'on', 'probeReading', scan.probe_reading, ... 
    'probeSensitivity', scan.probe_sensitivity, 'scanSensitivity', scan.scan_sensitivity, ... 
    'centerDataToPeak', 'on', 'micrometerUnits', micrometer_units);
[~, min_position] = min(scan.norm_trans);
zero_delay_shift_fs = scan.delay_fs(min_position);
scan.delay_fs = scan.delay_fs - zero_delay_shift_fs;

scan.norm_trans = 1 + scan.norm_trans;

end

