function [x_new,run_stat] = backtracking_line_search(runmod,x,delta_fx,delta_x,cf_params,alpha,beta,t_init)

if ~isvector(x),error('x has to be a vector.');end%if
if ~isvector(delta_fx),error('delta_fx has to be a vector.');end%if
if ~isvector(delta_x),error('delta_x has to be a vector');end %if
if size(x,2) ~= 1,x = x.';end%if
if size(delta_fx,2) ~= 1,delta_fx = delta_fx.';end%if
if size(delta_x,2) ~= 1,delta_x = delta_x.';end%if

if (length(x) ~= length(delta_fx)) | (length(x) ~= length(delta_x)) | (length(delta_fx) ~= length(delta_x))
    error('Lengths of the corresponding vectors have to be the same!');
end%if

if nargin < 6
    alpha = 0.1;
end%if
if nargin < 7
    beta = 0.7;
end%if
if isempty(alpha)
    alpha = 0.1;
end%if
if isempty(beta)
    beta = 0.7;
end%if

if nargin < 8
    t = 0.005/(norm(alpha*delta_fx.'*delta_x));
else
    t = t_init/(norm(alpha*delta_fx.'*delta_x));
end%if

eval_indicator = (delta_fx ~= 0);

f_x = evaluate_cost_function(x,cf_params,runmod);
%f_x = evaluate_partial_cost_function(x,eval_indicator,cf_params,runmod);
%f_x
number_of_iteration = 0;
end_of_iteration = 0;
number_of_inf = 0;

while end_of_iteration ~= 1
    number_of_iteration = number_of_iteration+1;
    
    [A,B,C,Q,sigma] = split_parameter_vector(x+t*delta_x);
    if any(Q ~= sort(Q))
        % if the monotonity is violated we return immediately
        t = t*beta;
        continue;
    end%if
    
    if sigma < 0
        t = t*beta;
        continue;
    end
     
    f_x_t_delta_x = evaluate_cost_function(x+t*delta_x,cf_params,runmod);  
    %f_x_t_delta_x = evaluate_partial_cost_function(x+t*delta_x,eval_indicator,cf_params,runmod);  
        %f_x_t_delta_x    
    if ~isreal(f_x_t_delta_x)
        disp('Not real..');
    end
    

    if isinf(f_x_t_delta_x)
        number_of_inf = number_of_inf+1;
    end%if
    if f_x_t_delta_x > (f_x + alpha*t*delta_fx.'*delta_x)
        t = t*beta;
    else
        end_of_iteration = 1;
    end%if
        
end%while end_of_iteration ~= 1

%disp(sprintf('%f',f_x_t_delta_x));
%disp(sprintf('number_of_iteration = %d, number_of_inf = %d',number_of_iteration,number_of_inf));

%We store the number of iterations to control the t_init.
run_stat.number_of_inf = number_of_inf;
run_stat.number_of_iteration = number_of_iteration;

x_new = x+t*delta_x;
