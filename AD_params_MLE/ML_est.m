function [A,B,C,tl,sigma] = ML_est(meas_data,runmod)

% This function calculates the maximum likelihood estimation.
%
% Input arguments:
%       meas_data : type of class mdata. Describes the measurement. For
%           example tha maximal analog transition level etc.
%       runmod : run modifier
%
% Output arguments:
%       A,B,C : result of the estimation w.r.t. parameters of the sine
%           wave. The assumed signal is
%           y(t) = A*sin(2*pi*freq*t)+B*cos(2*pi*freq*t)+C
%       tl : the transition levels
%       sigma : deviation of the additive Gaussian noise

if nargin < 2
    %in this case we have no runmod 
    runmod.verbose = 3;
    runmod.eps_limit = 0.01;
end%if;

Y_index = get_measured_data(meas_data);%get the indexes
freq = get_sine_freq(meas_data);
if isempty(freq)
    error('Four parameter fitting has not been implemented yet');
end%if
Ts = get_Ts(meas_data);

if runmod.verbose > 1
    disp('-----------------------------------------');
    disp('Estimation algorithm started.');
end%if

M = length(Y_index); %number of samples

%inital estimate of the transition levels
tl_ideal = get_tl_of_an_ideal_quantizer(meas_data);%these are the ideal transition levels, it is important that this step is the first



%Y_index = Y_index-min(Y_index);
meas_data = corrigate_minium_value_of_samples(meas_data);
[Ad,Bd,Cd] = init_sin_est(meas_data,runmod);%estimation in the digital domain
%during sine fitting the transitino levels are not used

A = d2c_amp(meas_data,Ad);
B = d2c_amp(meas_data,Bd);
C = d2c(meas_data,Cd);

%It is very important to remove missing code.
% To do this the indexes are modified and the transition levels are also
% corrected. tl will replace tl_ideal. It can be done because after this
% every parameter is in the analog domain.
Y_index_save = Y_index;
[Y_index,tl] = remove_missing_code(Y_index_save,tl_ideal);
meas_data = update_meas_data(meas_data,Y_index);

%Calculate how many different indecies are in Y_index
N = length(tl);

sigma = init_sigma_est(tl,meas_data,A,B,C,runmod);

P = construct_parameter_vector(A,B,C,tl,sigma,runmod);

cf_params.Ts = Ts;
cf_params.Y_index = Y_index;
cf_params.freq = freq;
cf_params.meas_data = meas_data;

tl_grad = zeros(N,1);%init

eps_limit = runmod.eps_limit;
end_of_iteration = 0;

t_init_tl = 0.05*ones(size(tl));

%input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
input_signal = reconstruct_input_signal(meas_data,A,B,C);

prob = calculate_probabilities(tl,sigma,Y_index,input_signal,runmod);% size = M x 1
if any(~(prob >0))
    disp('**Problem: Invalid initial state.');
end

while(~end_of_iteration)
    P_prev = P;

    %---------------------------------------------------------------------
    %calculate the gradient wrt. A,B,C -- the input signal
    [A_grad,B_grad,C_grad] = gradient_wrt_ABC(A,B,C,tl,sigma,meas_data,runmod);
    P_grad_ABC = ABC_part_of_gradient_vector(A_grad,B_grad,C_grad,length(P));

    %calculate the backtracking linesearch
    if runmod.verbose >= 1, disp('Backtracking w.r.t ABC');end;
    P_temp = P;
    P = backtracking_line_search(runmod,P,P_grad_ABC,-P_grad_ABC,cf_params,0.005);

    if runmod.verbose >= 1, disp(sprintf('delta ABC = %f',num2lsb(meas_data,norm(P-P_temp))));end;
    [A,B,C,tl,sigma] = split_parameter_vector(P);

    P_prev_temp = P;

    %---------------------------------------------------------------------    

    [P,t_init_tl] = gradest_tl(P,Y_index,Ts,freq,meas_data,t_init_tl,runmod);

    %---------------------------------------------------------------------

    if runmod.verbose >= 1, disp('Backtracking w.r.t sigma');end;
    
    sigma_grad = gradient_wrt_sigma(A,B,C,tl,sigma,meas_data,runmod);
    grad_P_sigma = sigma_part_of_gradient_vector(sigma_grad,length(P));

    P_temp = P;
    P = backtracking_line_search(runmod,P,grad_P_sigma,-grad_P_sigma,cf_params,0.6,0.7,10^4);
    disp(sprintf('delta sigma = %f',num2lsb(meas_data,norm(P-P_temp))));

    %---------------------------------------------------------------------

    %if norm(P-P_prev) < eps_limit
    if num2lsb(meas_data,max(abs(P-P_prev))) < eps_limit
        end_of_iteration = 1;%stop the iteration
    end%if

end%while(~end_of_iteration)

% if strcmp(runmod.tl_alg,'groupped')
%     rm_temp = runmod;
%     rm_temp.tl_alg = 'full';
%     end_of_iteration = 0;
%     while(~end_of_iteration)
%         [P,t_init_tl] = gradest_tl(P,Y_index,Ts,freq,meas_data,t_init_tl,rm_temp);
%         if num2lsb(meas_data,max(abs(P-P_prev))) < eps_limit
%             end_of_iteration = 1;%stop the iteration
%         end%if        
%     end%while(~end_of_iteration)
% end%if

[A,B,C,tl,sigma] = split_parameter_vector(P);

