% This script establishes an environment to run the Maximum Likelihood Estimation script
% under development.  

close all;
clear all;
clc;

diary ML_est.log
diary on;

Vmax =  2;
Vmin = -2;

rm.verbose = 2; %run modifier
%rm.verbose = 3; %run modifier
rm.eps_limit = 0.05;
rm.Tl_alg = 'groupped';
%rm.tl_alg = 'full';


MC_run = 10; %Number of Monte Carlo runs
N=8;          %Number of bits

%Initial values for the input sine wave and the noise
freq = 1;
dc_level = 0;
amplitude = 2.021;
%amplitude = 1.613;      %This is for calculating the midcode
Ts = 1/10e3;
N_samples = 8*80000;
sigma = .018;


rm.noise_model = 'Gauss';
MC_results = zeros(MC_run,5);



for k=1:MC_run

    disp(sprintf('Monte Carlo Iteration #: %d ',k));

    meas_data = mdata;
    meas_data = set_V_max(meas_data,Vmax);
    meas_data = set_V_min(meas_data,Vmin);
    meas_data = set_Nbit(meas_data,N);
    meas_data = set_Ts(meas_data,Ts);
    meas_data = set_sine_freq(meas_data,freq);

    Tl_ideal = get_Tl_of_an_ideal_quantizer(meas_data);

    tic;            %start timer
    toc_temp = toc;
    [bb,cc] = gen_samples(Tl_ideal,sigma,amplitude,freq,dc_level,Ts,N_samples);

    NN = length(bb);
    tt_bb = 0:NN-1;
    tt_bb = tt_bb.';
    bb2 = [bb(1:end/4); bb(3*end/4:end)];
    tt_bb2 = [tt_bb(1:end/4); tt_bb(3*end/4:end)];
    meas_data = set_meas_data(meas_data,bb,tt_bb);

    %%%%%%%%%%%%%%%%%%%%%%%% Plot/Debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rm.verbose > 2
        figure; plot(tt_bb,bb); grid; hold; plot(20*cc,'r');
    end%if
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    code_density_INL_DNL(meas_data,rm);
    toc_save = toc-toc_temp;

%     MC_results(k,:) = [A_est,B_est,C_est,sigma_est,toc_save];
    close all;
end