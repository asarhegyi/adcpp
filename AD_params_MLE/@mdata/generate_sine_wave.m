function [S,C] = generate_sine_wave(obj,freq)

%sfreq = get_sine_freq(meas_data);
if isempty(obj.sine_freq) & (nargin < 2)
    error('Reconstructing the input signal requires the frequency of the sine wave');
end%if

if nargin < 2
    freq = obj.sine_freq;
end

tt = obj.measure_time;
if isempty(tt)
    M = length(obj.measured_data);
    tt = 0:M-1;
    tt = tt.';
end%if

tt = tt*obj.Ts;

S = sin(2*pi*freq*tt);
C = cos(2*pi*freq*tt);
