% This script runs the Monte Carlo Simulation for the ML estimation
% the rm.noise_model modifier selects the noise model discussed in Chapter
% 6 of the thesis.

close all;
clear all;
clc;

Vmax =  2;
Vmin = -2;

rm.verbose = 2; %run modifier
%rm.verbose = 3; %run modifier
rm.eps_limit = 0.05;
rm.tl_alg = 'groupped';
%rm.tl_alg = 'full';

%Select noise model
%Noise_models = 'Gauss';    %Gaussian noise model
%Noise_models = 'Laplace';  %Laplacian noise model
%Noise_models = 'Triangle'; %Triangular noise model

%Noise_models_lst = ['Gauss   ';'Laplace ';'Triangle'];
Noise_models_lst = ['Gauss   ';'Laplace '];
Noise_models = cellstr(Noise_models_lst);

M = 100;    %Number of Monte Carlo runs
N=8;        %Number of bits

%Initial values for the input sine wave and the noise
freq = 1;
dc_level = 0;
amplitude = 2.021;
Ts = 1/10e3;
N_samples = 8*8000;
sigma = .018;


for model_cnt=1:length(Noise_models) 

    rm.noise_model = char(Noise_models(model_cnt));
    MC_results = zeros(M,3+2^N-1+1+1);

    for k=1:M

        disp(sprintf('Iteration #: %d ',k));

        meas_data = mdata;
        meas_data = set_V_max(meas_data,Vmax);
        meas_data = set_V_min(meas_data,Vmin);
        meas_data = set_Nbit(meas_data,N);
        meas_data = set_Ts(meas_data,Ts);
        meas_data = set_sine_freq(meas_data,freq);

        tl_ideal = get_tl_of_an_ideal_quantizer(meas_data);

        tic;            %start timer
        toc_temp = toc;
        [bb,cc] = gen_samples(tl_ideal,sigma,amplitude,freq,dc_level,Ts,N_samples);

        NN = length(bb);
        tt_bb = 0:NN-1;
        tt_bb = tt_bb.';
        bb2 = [bb(1:end/4); bb(3*end/4:end)];
        tt_bb2 = [tt_bb(1:end/4); tt_bb(3*end/4:end)];
        meas_data = set_meas_data(meas_data,bb2,tt_bb2);

        [A_est,B_est,C_est,Tl_est,sigma_est] = ML_est(meas_data,rm);
        toc_save = toc-toc_temp;

        MC_results(k,:) = [A_est,B_est,C_est,Tl_est',sigma_est,toc_save];

    end

    disp(sprintf('MC_results_%s.txt',rm.noise_model));
    FILEO = sprintf('MC_results_%s.txt',rm.noise_model);
    save(FILEO,'MC_results','-ASCII'); 

end