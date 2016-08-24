function tl = code2tl(a,Q)

% Input check

%Size check 
if ~ ( ( size(Q,1) == 1) | (size(Q,2) == 1))
    error('Q must a vector');
end%if

if length(Q) < 2
    error('Length of Q must be equal to or greater than 2');
end%if

% it must a sorted sequence
if Q ~= sort(Q)
    error('Q must be sorted');
end%if

if Q(1) ~= min(Q) 
    error('Q must be in ascending order');
end%if

% do we have enough information avaible?
if isempty(a.analog_min)     
    error('In the measurement description object the analog minimum value is missing.');
end%if

if isempty(a.analog_max)
    error('In the measurement description object the analog maximum value is missing.');
end

Q_min = min(Q);
Q_max = max(Q);

% if ~isempty(a.digital_min)
%     if Q_min ~= a.digital_min
%         warning('Digital minimum will be overwritten');
%     end
% end%if
% 
% if ~isempty(a.digital_max)
%     if Q_max ~= a.digital_max
%         warning('Digital maximum will be overwritten');
%     end
% end%if

a.digital_min = Q_min;
a.digital_max = Q_max;

if length(Q) < length(Q_min:1:Q_max)
    error('Missing code is not supported yet');
end


%tl = (Q-Q_min)/(Q_max-Q_min); %normalize
if length(Q) > 2
    tl = (0:length(Q)-2)/(length(Q)-2);
    tl = tl*(a.analog_max-a.analog_min) + a.analog_min;
else
    tl = mean([a.analog_min, a.analog_max]);
end%if

tl = tl.';

