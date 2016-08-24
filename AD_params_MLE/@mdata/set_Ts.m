function obj = set_Ts(obj,Ts_in)

if ~(Ts_in > 0)
    error('Sampling time must be positive');
end%if

obj.Ts = Ts_in;
