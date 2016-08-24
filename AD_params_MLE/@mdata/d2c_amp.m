function p=d2c_amp(obj,pd)

% if isempty(obj.analog_min)
%     error('Analog minimum level is missing')
% end
% 
% if isempty(obj.analog_max)
%     error('Analog maximum level is missing')
% end
% 
% if isempty(obj.digital_min)
%     error('Digital minimum level is missing')
% end
% 
% if isempty(obj.digital_max)
%     error('Digital maximum level is missing')
% end

% y = ax+b
% x:analog domain; y:digital domain
% x1 = analog_max, y1 = digital_max
% x2 = analog_min, y2 = digital_min

% a = (y1-y2)/(x1-x2)
% b = (y2x1-y1x2)/(x1-x2)

N = 2.^(obj.N_bit);
N = obj.N;

%a = (obj.digital_max-obj.digital_min)/(obj.analog_max-obj.analog_min);
a = N/(obj.V_max-obj.V_min);

p = pd/a;
