function data_obj = mdata()

% Constructor method of the object mdata.
% Usage:
%       data_obj = mdata()
% where
%       data_obj : the result of the constructor method.
%

data_obj.V_min = []; %analog value of the least digital output
data_obj.V_max = []; %analog value of the largest digital output

data_obj.N_bit = []; % number of bits

data_obj.N = []; % number of output code

%data_obj.digital_min = [];
%data_obj.digital_max = [];

data_obj.measured_data = [];
data_obj.measure_time = []; %if it is not empty then time is not in form [0:M-1]*Ts
data_obj.Ts = [];%sampling time
data_obj.sine_freq = []; %in the case 3 parameter sine fitting it must be given

data_obj = class(data_obj,'mdata');
