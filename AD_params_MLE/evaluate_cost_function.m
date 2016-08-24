function value = evaluate_cost_function(P,cf_params,runmod)

if ~isvector(P)
    error('P has to be a vector.');
end 

if size(P,2) ~= 1
    P = P.';
end%if

if length(P) < 5
    error('P is too short.');
end%if

[A,B,C,Q,sigma] = split_parameter_vector(P);

Ts = cf_params.Ts;
freq = cf_params.freq;
Y_index = cf_params.Y_index;
M = length(Y_index);

%input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
input_signal = reconstruct_input_signal(cf_params.meas_data,A,B,C);

% calculate the probalities
prob = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod);

if(prob<0)
    error('Probability is less than 0!');
end%if(prob<0)

% little trcik to avoid warnings
prob_zero = find(prob == 0.0);

if isempty(prob_zero)
    value = -sum(log(prob));    
else
    value = Inf;
    
    %debug
    %warning('Probability is corrected!');
    prob(prob_zero) = 10^-7;
    value = -sum(log(prob));
end%if



