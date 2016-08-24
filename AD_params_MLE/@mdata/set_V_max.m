function a = set_V_max(a,value)

if ~isempty(a.V_min)
    if ~(a.V_min < value)
        error('Maximum value must be greater than mininum value'); 
    end%if
end%if


a.V_max = value;
