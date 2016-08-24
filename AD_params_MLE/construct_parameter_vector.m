function P = construct_parameter_vector(A,B,C,Q,sigma,runmod)

if nargin < 6
    runmod.verbose = 3;
end%if;

P = zeros(length(Q)+4,1);

P(1) = A;
P(2) = B;
P(3) = C;
P(4:end-1) = Q;
P(end) = sigma;
