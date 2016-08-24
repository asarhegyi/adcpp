function [A,B,C,Q,sigma] = split_parameter_vector(P)

A = P(1);
B = P(2);
C = P(3);
Q = P(4:end-1);
sigma = P(end);
