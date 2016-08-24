function obj = set_meas_data(obj,data,time_data)

obj.measured_data = data;

if nargin > 2
    obj.measure_time = time_data;
else
    obj.measure_time = [];
end%if

