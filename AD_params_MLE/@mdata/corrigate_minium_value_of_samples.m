function obj = corrigate_minium_value_of_samples(obj)


if ~isempty(obj.measured_data)
    obj.measured_data = obj.measured_data-min(obj.measured_data);
end%if
