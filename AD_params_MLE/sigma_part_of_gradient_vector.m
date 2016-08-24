function [grad_P_sigma] = sigma_part_of_gradient_vector(grad_sigma,length_of_P)

grad_P_sigma = zeros(length_of_P,1);

grad_P_sigma(end) = grad_sigma;
