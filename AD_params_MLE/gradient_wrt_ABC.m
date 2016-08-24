function [A_grad,B_grad,C_grad] = gradient_wrt_ABC(A,B,C,Q,sigma,meas_data,runmod)

% calculates the gradient wrt the parameters of the sine wave

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

numerator = zeros(M,1);

if strcmp(runmod.noise_model,'Gauss')
    const_outer = 1/( sqrt(pi*2)*(sigma));
    const_inner = -1/(2*(sigma^2));

    numerator(Y_index_others) = const_outer *( exp(const_inner*( temp_others_k.^2)) - exp(const_inner*( temp_others_k_minus_1.^2)));
    numerator(Y_index_0) = const_outer * ( exp(const_inner*( temp_0.^2)) );
    numerator(Y_index_N) = -const_outer * ( exp(const_inner*( temp_N.^2)) );

elseif strcmp(runmod.noise_model,'Laplace')
    const_outer = 1/( 2*sigma );
    const_inner = -1/sigma;

    numerator(Y_index_others) = const_outer *( exp(const_inner*abs(temp_others_k)) - exp(const_inner*abs(temp_others_k_minus_1)) );
    numerator(Y_index_0) = const_outer * ( exp(const_inner*abs(temp_0)) );
    numerator(Y_index_N) = -const_outer * ( exp(const_inner*abs(temp_N)) ); 
end%if

denominator = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod);% size = M x 1 

d_0 = find(denominator == 0);
denominator(d_0) = 0.00001;
if ~isempty(d_0)
    warning('Probability is corrected [gradient_wrt_ABC]');
end

% tt = 0:M-1;
% tt = tt.';
% tt = tt*Ts;
[S,C] = generate_sine_wave(meas_data); %S: sin(2*pi*freq*tt), C: cos(2*pi*freq*tt), 

% - because of minimization, ./denominator because of log likelihood cost function
% A_grad = -sum(numerator./denominator.*(-sin(2*pi*freq*tt)));% second -  because temp = (Q()-input_signal)
% B_grad = -sum(numerator./denominator.*(-cos(2*pi*freq*tt)));
% C_grad = -sum(numerator./denominator*(-1));
% A_grad = -mean(numerator./denominator.*(-sin(2*pi*freq*tt)));% second -  because temp = (Q()-input_signal)
% B_grad = -mean(numerator./denominator.*(-cos(2*pi*freq*tt)));
% C_grad = -mean(numerator./denominator*(-1));

A_grad = -mean(numerator./denominator.*(-S));% second -  because temp = (Q()-input_signal)
B_grad = -mean(numerator./denominator.*(-C));
C_grad = -mean(numerator./denominator*(-1));
