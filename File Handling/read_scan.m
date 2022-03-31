function [x_data_mm_or_fs, signal] = read_scan(varargin)

p = inputParser;
p.addParameter('filename', '');
p.addParameter('channel', 1);
p.addParameter('micrometerUnits', 'mm');
p.addParameter('convertToDelay', 'on');
p.addParameter('plot', 'on');
p.addParameter('centerDataToPeak', 'on');
p.addParameter('normalizeToProbe', 'off');
p.addParameter('probeSensitivity', 0);
p.addParameter('probeReading', 0);
p.addParameter('scanSensitivity', 0);

parse(p, varargin{:})
filename = p.Results.filename;
plot_results_flag = p.Results.plot;
micrometer_units = p.Results.micrometerUnits;
convert_to_delay = p.Results.convertToDelay;
center_data_to_peak = p.Results.centerDataToPeak;
probe_sensitivity = p.Results.probeSensitivity;
probe_reading = p.Results.probeReading;
scan_sensitivity = p.Results.scanSensitivity;
normalize_to_probe = p.Results.normalizeToProbe;

if ~isempty(filename)
    filename = file_selector('filename', filename);
else
    filename = file_selector();
end

[stage_position, signal] = get_scan_data(filename);

if strcmpi(micrometer_units, 'in')
    stage_position = stage_position*2.54;
end

if strcmpi(center_data_to_peak, 'on')
    [stage_position, signal] = center_data_to_peak_fcn(stage_position, signal);
end

if strcmpi(convert_to_delay, 'on')
    x_data_mm_or_fs = 2*stage_position/1000/2.998e8*1e15;
else
    x_data_mm_or_fs = stage_position;
end

if strcmpi(normalize_to_probe, 'on')
    if probe_sensitivity == 0 || scan_sensitivity == 0 || probe_reading == 0
        error('If normalizeToProbe option is selected, you must provide the probe sensitivity, probe reading and scan sensitivity')
    end
    scaled_probe_reading = probe_sensitivity*probe_reading;
    signal = signal*scan_sensitivity/scaled_probe_reading;
end
if strcmpi(plot_results_flag, 'on')
    plot(x_data_mm_or_fs, signal, 'r')
    
    if strcmpi(convert_to_delay, 'on')
        xlabel('delay(fs)');
    else
        xlabel('stage_position (mm)');
    end
        ylabel('signal (arb)');
end
