function [Y_index] = quantize_samples(Q,signal_in)

if( size(Q,1) ~= 1)
    Q = Q.';
end%if

% slack variables to construct the quatized signal
%zz2 = zz*ones(1,length(Q(1:end)));
%qq2 = ones(size(zz))*Q(1:end);

%yy = sign_modified(zz2-qq2);%% yy the indexes of the decision levels
%yy = sum(yy,2);

%Y_index = yy;

Y_index = quant(Q,signal_in);

function y = sign_modified(x)
y = sign(x)+(x==0);
y = (y+1)*.5;
