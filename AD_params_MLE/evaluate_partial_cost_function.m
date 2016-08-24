function value = evaluate_partial_cost_function(P,P_eval_indicator,cf_params,runmod)

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

%check which parameter will change (in the backtracking algorithm)
[ind_A,ind_B,ind_C,ind_Q,ind_sigma] = split_parameter_vector(P_eval_indicator);

Ts = cf_params.Ts;
freq = cf_params.freq;
Y_index = cf_params.Y_index;
M = length(Y_index);

input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);


% If A,B,C or sigma are in the set of parameters which are modified during backtracking, then we cannot short by partial evaluate the cost function.
% First, check if at least one of the parameters above mentioned is zero or not.
int_total_evaluation = isempty(find(ind_A ~= 0))*isempty(find(ind_B ~= 0))*isempty(find(ind_C ~= 0))*isempty(find(ind_sigma ~= 0));


% calculate the probalities
if int_total_evaluation == 0
    prob = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod);
else
     Q_temp_1 = [(ind_Q(1:end-1) ~= 0);0];
     Q_temp_2 = [0;(ind_Q(1:end-1) ~= 0)];
     Q_temp = Q_temp_1+Q_temp_2;
     if(ind_Q(end) ~= 0)
         Q_temp(end) = 1;
     end%if     
     
     q = Q_temp.*[1:length(Q_temp)].';
     %QQ = Q_temp *ones(1,M);
     QQ = q*ones(1,M);
     YY = ones(length(Q),1)*(Y_index.');
     %Y_index_partial = find(Q_temp > 0);
     
     if(ind_Q(end) ~= 0)
         QQ_tail = (q(end)+1)*ones(size(Q_temp))*ones(1,M);
         Y_index_partial = find(sum( (QQ == YY) + (QQ_tail == YY)));
     else
         Y_index_partial = find(sum(QQ == YY));
     end
    prob = calculate_probabilities(Q,sigma,Y_index(Y_index_partial),input_signal(Y_index_partial),runmod);
end %if

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
    warning('Probability is corrected!');
    prob(prob_zero) = 10^-7;
    value = -sum(log(prob));
end%if



