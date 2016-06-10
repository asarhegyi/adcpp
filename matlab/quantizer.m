%% MATLAB
function out = quantizer(args)


% This script intents to verify the Cpp quantizer. It generates
% the code transition levels for an ideal quantizer. Vmax, Vmin, and
% qbits parameters have to match with the ones in gen_stimulus.cpp.
% This has to be done manually. First gen_stimulus.cpp has to compiled
% and run, which generates among others noisyData.dat. This script grabs
%  noisyData.dat and quantize it with the matlab quantizer. 
% Both the code transition levels and the results of the quantization are
% saved to .dat files which can be compared with the Cpp results.
%
%       out = quantizer(args)
%
%       Input arguments:
%         args.arg1: number of quantizer bits
%
%       This program is public domain. It is available through
%       https://github.com/asarhegyi/

% $Id:  $
% Copyright (c) 2015-2016 by Attila Sarhegyi
% All rights reserved.


qbits=args.arg1;

Vmax =  2.15;
Vmin = -2.15;

meas_data = mdata;
meas_data = set_V_max(meas_data,Vmax);
meas_data = set_V_min(meas_data,Vmin);
meas_data = set_Nbit(meas_data,qbits);

load noisyData.dat

tl_ideal = get_Tl_of_an_ideal_quantizer(meas_data);
Y_index = quantize_samples(tl_ideal,noisyData);

% save the code transition levels and compare with Cpp
fid = fopen('Tl.mat.dat', 'w');
fprintf(fid, '%.14f\n', tl_ideal);
fclose(fid);

% save the code transition levels and compare with Cpp
fidq = fopen('quantizedData.mat.dat', 'w');
fprintf(fidq, '%d\n', Y_index);
fclose(fidq);
