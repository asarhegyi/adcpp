%% MATLAB
function ret = wrapper_sfit(args)


% This wrapper provides and interface between the python
% testbench (python-matlab-bridge) and the MATLAB sine-fit
% scripts (sfit3.m, sfit4imp.m) for the Monte Carlo simulation.
%
%       ret = wrapper_sfit(args)
%        
%       Input arguments:
%         args.arg1: the input file name with a vector of samples
%         args.arg2: number of parameters to estimate (3 or 4)
%         args.arg3: sample rate
%         args.arg4: fundamental frequency for sfit3 and
%                starting frequency for sfit4imp
%
%       This program is public domain. It is available through
%       https://github.com/asarhegyi/

% $Id:  $
% Copyright (c) 2015-2016 by Attila Sarhegyi
% All rights reserved.


    %start timer
    tic;

    data = load(args.arg1);
    numParams = args.arg2;
    sample_rate = args.arg3;
    input_fr = args.arg4;

    time = [0:1:length(data)-1]';
    time = time*1/sample_rate;

    if numParams == 3
        [X, Rn] = sfit3(data, time, sample_rate, input_fr);
    elseif numParams == 4
        %[X, Rn] = sfit4(data, time, sample_rate);
        [X, Rn] = sfit4imp(data, '', sample_rate, '', 'Yes');
    else
        error('Invalid numParams argument');
    end

    %stop timer
    toc_save = toc;

    ret.exec_time = sprintf('%.6f', toc_save);
    %disp(ret.exec_time);

    ret.sample_rate  = sprintf('%.20f', sample_rate);
    %disp(ret.sample_rate);

    ret.A  = sprintf('%.20f', X.A);
    ret.f  = sprintf('%.20f', X.f);
    ret.phi= sprintf('%.20f', X.phi*(pi/180));
    ret.dc = sprintf('%.20f', X.DC);
    ret.erms = sprintf('%.20f', X.erms);

    %disp(ret.A);
    %disp(ret.f);
    %disp(ret.phi);
    %disp(ret.dc);
    %disp(ret.erms);

end
