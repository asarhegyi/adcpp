function check_indecies(Y_index)

% Function check_indecies(Y_index) checks if the given set of measured
% codes consitutue a valid set.  
%   - It must contain integer values
%   - The least value must be 1
% The function throws error if at least one of the items from the previous
% list does not come true (see the detailed error message). 

%  check int property
if Y_index ~= floor(Y_index)
    error('Output indecies of the AD must be integer numbers');
end%if

if min(Y_index) ~= 0
    error('Minimum index must be 0');
end%if
