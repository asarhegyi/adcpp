function [L,grad] = Likelihood_wrt_ABCN(x,Tl,meas_data,runmod)

N = length(Tl);         %number of quantazation levels
Y_index = get_measured_data(meas_data);
M = length(Y_index);                        %number of samples

Y_index_0 = find(Y_index == 0);
Y_index_N = find(Y_index == N);
Y_index_others = find( ((Y_index ~= 0) & (Y_index ~= N) ));

%Using the given parameteres we reconstruct the input signal
if strcmp(runmod.sine,'AphiDC')
    input_signal = reconstruct_input_signal2(meas_data,x(1),x(2),x(3));
else
    input_signal = reconstruct_input_signal(meas_data,x(1),x(2),x(3));
end
sigma = x(4);

%%%%%%%%%%%%%%%%%%%%%%% calculate Likelihood function %%%%%%%%%%%%%%%%%%%%%
prob = zeros(M,1);
    
if ~isempty(Y_index_others)
    prob(Y_index_others) = .5*(erf( ( Tl(Y_index(Y_index_others)+1) - input_signal(Y_index_others))/sigma/sqrt(2)) ...
        - erf(( Tl(Y_index(Y_index_others)) - input_signal(Y_index_others))/sigma/sqrt(2)));
end%if

if ~isempty(Y_index_0)
    prob(Y_index_0) = 0.5*(1+erf( ( Tl(1) - input_signal(Y_index_0))/sigma/sqrt(2)) );
end%if

if ~isempty(Y_index_N)
    prob(Y_index_N) = 0.5*(1-erf( ( Tl(N) - input_signal(Y_index_N))/sigma/sqrt(2)));
end%if

if any(~(prob >0))
    disp('**Problem: Zero Probability.');
    d_0 = logical(prob == 0);
    prob(d_0) = eps;
end

logprob = log(prob);
L = -1*sum(logprob,1);
    

%%%%%%%%%%%%%%%%%%%%%%% calculate gradients wrt ABC %%%%%%%%%%%%%%%%%%%%%%%

temp_0 = (Tl(1)-input_signal(Y_index_0));
temp_N = (Tl(N)-input_signal(Y_index_N));

temp_others_k = (Tl(Y_index(Y_index_others)+1)-input_signal(Y_index_others));
temp_others_k_minus_1 = (Tl(Y_index(Y_index_others))-input_signal(Y_index_others));

numerator = zeros(M,1);
grad_sigma = zeros(M,1);

const_outer = 1/( sqrt(pi*2)*(sigma));
const_inner = -1/(2*(sigma^2));

numerator(Y_index_others) = const_outer *( exp(const_inner*( temp_others_k.^2)) - ...
    exp(const_inner*( temp_others_k_minus_1.^2)));
numerator(Y_index_0) = const_outer * ( exp(const_inner*( temp_0.^2)) );
numerator(Y_index_N) = -const_outer * ( exp(const_inner*( temp_N.^2)) );

const_outer = 1/( sqrt(2*pi)*(sigma^2));    %redefine multipier for grad_sigma

%calculate dP/d_sigma
grad_sigma(Y_index_0) = const_outer*exp(const_inner*(temp_0.^2)) .* (-temp_0);
grad_sigma(Y_index_N) = const_outer*exp(const_inner*(temp_N.^2)) .* (temp_N);
grad_sigma(Y_index_others) = const_outer*exp(const_inner*(temp_others_k.^2)) .* (-temp_others_k) + ...
    const_outer*exp(const_inner*(temp_others_k_minus_1.^2)) .* (temp_others_k_minus_1);

denominator = prob;

% d_0 = find(denominator == 0);
% denominator(d_0) = 0.00001;
% if ~isempty(d_0)
% %    warning('Probability is corrected [Likelihood_function_wrt_ABC]');
%     disp('**Probability is corrected.');
% end

[S,C] = generate_sine_wave(meas_data);      %S: sin(2*pi*freq*tt), C: cos(2*pi*freq*tt), 

if strcmp(runmod.sine,'AphiDC')
    A_grad = sum(numerator./denominator.*(-C));
    B_grad = sum(numerator./denominator.*(x(1)*S));
    C_grad = sum(numerator./denominator*(-1));
else
    A_grad = sum(numerator./denominator.*(-S));
    B_grad = sum(numerator./denominator.*(-C));
    C_grad = sum(numerator./denominator*(-1));
end

sigma_gradient = sum(grad_sigma./denominator);

if nargout > 1
grad(1)=A_grad;
grad(2)=B_grad;
grad(3)=C_grad;
grad(4)=sigma_gradient;
end