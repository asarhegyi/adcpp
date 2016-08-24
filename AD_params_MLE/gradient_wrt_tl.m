function [Q_grad] = gradient_wrt_tl(dQ_index,A,B,C,Q,sigma,meas_data,runmod)
    
Y_index = get_measured_data(meas_data);    
% calculates the gradient wrt a quantization level

N = length(Q);%N = number of quantazation levels
% => range(Y_index) = 0,...,N

M = length(Y_index);%number of samples

Y_index_k = find(Y_index == dQ_index-1);
Y_index_k_plus_1 = find(Y_index == dQ_index);

%input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
input_signal = reconstruct_input_signal(meas_data,A,B,C);

temp_k = (Q(dQ_index)-input_signal(Y_index_k));
temp_k_plus_1 = (Q(dQ_index)-input_signal(Y_index_k_plus_1));

numerator = zeros(M,1);

if strcmp(runmod.noise_model,'Gauss')
    const_outer = 1/( sqrt(pi*2)*(sigma));
    const_inner = -1/(2*(sigma^2));

    numerator(Y_index_k) = const_outer *( exp(const_inner*( temp_k.^2)));
    numerator(Y_index_k_plus_1) = const_outer * (- exp(const_inner*( temp_k_plus_1.^2)));

elseif strcmp(runmod.noise_model,'Laplace')
    % This is fitted to the Laplace noise model
    const_outer = 1/( 2*sigma );
    const_inner = -1/sigma;

    numerator(Y_index_k) = const_outer *( exp(const_inner*abs( temp_k )));
    numerator(Y_index_k_plus_1) = const_outer * (- exp(const_inner*abs( temp_k_plus_1 ))); 
end%if

denominator = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod);% size = M x 1

d_0 = find(denominator == 0);
denominator(d_0) = 0.0000001;
if ~isempty(d_0)
    warning('Probability is corrected [gradient_wrt_tl]');
end

Q_grad = -mean(numerator./denominator);% - because of minimization, ./denominator because of log likelihood cost function
