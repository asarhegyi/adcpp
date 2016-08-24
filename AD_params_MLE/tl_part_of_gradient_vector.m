function [grad_P_tl] = tl_part_of_gradient_vector(grad_tl,tl_index,length_of_P)

grad_P_tl = zeros(length_of_P,1);

grad_P_tl(3+tl_index) = grad_tl;
