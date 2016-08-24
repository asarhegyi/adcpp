function [grad_P_ABC] = ABC_part_of_gradient_vector(A_grad,B_grad,C_grad,length_of_P)

grad_P_ABC = zeros(length_of_P,1);

grad_P_ABC(1) = A_grad;
grad_P_ABC(2) = B_grad;
grad_P_ABC(3) = C_grad;

