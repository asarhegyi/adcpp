function input_signal = reconstruct_input_signal(obj,A,phi,DC,freq)

%sfreq = get_sine_freq(meas_data);
if isempty(obj.sine_freq) & (nargin < 5)
    error('Reconstructing the input signal requires the frequency of the sine wave');
end%if

if nargin < 5
    freq = obj.sine_freq;
end

tt = obj.measure_time;
if isempty(tt)
    M = length(obj.measured_data);
    tt = 0:M-1;
    tt = tt.';
end%if

tt = tt*obj.Ts;
%input_signal = A*sin(2*pi*freq*tt)+B*cos(2*pi*freq*tt)+C;
input_signal = A*cos(2*pi*freq*tt+phi)+DC;



