function sigma_gradient = gradient_wrt_sigma(A,B,C,Q,sigma,meas_data,runmod)

% Here we assume that the cost function should be minimized (not
% maximized).

Y_index = get_measured_data(meas_data);

N = length(Q);%N = number of quantazation levels
% => range(Y_index) = 1,...,N+1

M = length(Y_index);%number of samples

Y_index_0 = find(Y_index == 0);
Y_index_N = find(Y_index == N);
Y_index_others = find( ((Y_index ~= 0) & (Y_index ~= N) ));

%input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
input_signal = reconstruct_input_signal(meas_data,A,B,C);

temp_0 = (Q(1)-input_signal(Y_index_0));
temp_N = (Q(N)-input_signal(Y_index_N));

temp_others_k = (Q(Y_index(Y_index_others)+1)-input_signal(Y_index_others));
temp_others_k_minus_1 = (Q(Y_index(Y_index_others))-input_signal(Y_index_others));

grad_sigma = zeros(M,1);

if strcmp(runmod.noise_model,'Gauss')
    %we define two constants in order to simplify the computation
    const_outer = 1/( sqrt(2*pi)*(sigma^2));
    const_inner = -1/(2*(sigma^2));

    %calculate dP/d sigma
    grad_sigma(Y_index_0) = const_outer*exp(const_inner*(temp_0.^2)) .* (-temp_0);
    grad_sigma(Y_index_N) = const_outer*exp(const_inner*(temp_N.^2)) .* (temp_N);
    grad_sigma(Y_index_others) = const_outer*exp(const_inner*(temp_others_k.^2)) .* (-temp_others_k) + ...
       const_outer*exp(const_inner*(temp_others_k_minus_1.^2)) .* (temp_others_k_minus_1);

elseif strcmp(runmod.noise_model,'Laplace')
    % This is fitted to the Laplace noise model
    const_outer = 1/( 2*sigma );
    const_inner = -1/sigma;

    grad_sigma(Y_index_0) = const_outer*exp(const_inner*abs( temp_0 )) .* (-temp_0);
    grad_sigma(Y_index_N) = const_outer*exp(const_inner*abs( temp_N )) .* (temp_N);
    grad_sigma(Y_index_others) = const_outer*exp(const_inner*abs( temp_others_k )) .* (-temp_others_k) + ...
        const_outer*exp(const_inner*abs( temp_others_k_minus_1 )) .* (temp_others_k_minus_1); 
end%if

prob = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod);% size = M x 1

sigma_gradient = -sum(grad_sigma./prob);% - because of minimization, ./prob because of log likelihood cost function
