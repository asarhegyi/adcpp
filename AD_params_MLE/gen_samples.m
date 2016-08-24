function [Y_index,original_signal,original_signal_with_noise] = gen_samples(Q,sigma,A,freq,dc,Ts,N)

% Q : quantization levels
% sigma^2 : additive noise variance
% A : amplitude of the sine
% freq : frequency of the sine
% dc : dc level of the sine
% Ts : sampling time    \
% N : number of the generated samples

%we construct the noisy signal
tt = 0:N-1;
tt = tt*Ts;
zz = A*sin(2*pi*freq*tt+pi*1.2)+dc;

zz_noise_free = zz;

zz = zz.';
zz = zz+sigma*randn(size(zz));

Y_index = quantize_samples(Q,zz);

original_signal = zz_noise_free;

original_signal_with_noise = zz;
