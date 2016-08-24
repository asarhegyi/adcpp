function obj = update_meas_data(obj,data)

if size(data) ~= size(obj.measured_data)
    error('Invalid measurement update. Sizes must be correct.');
end

obj.measured_data = data;

