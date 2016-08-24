function a = set_V_min(a,value)

if ~isempty(a.V_max)
    if ~(a.V_max > value)
        error('Maximum value must be greater than mininum value'); 
    end%if
end%if


a.V_min = value;
