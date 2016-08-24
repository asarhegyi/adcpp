function obj = set_sine_freq(obj,freq)

if ~(freq > 0)
    error('Frequency must be positive');
end%if

obj.sine_freq = freq;
