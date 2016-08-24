function prob = calculate_probabilities(Q,sigma,Y_index,input_signal,runmod)

N = length(Q);%N = number of quantazation levels
% => range(Y_index) = 0,...,N

M = length(Y_index);%number of samples

Y_index_0 = find(Y_index == 0);
Y_index_N = find(Y_index == N);
Y_index_others = find( ((Y_index ~= 0) & (Y_index ~= N) ));

%tt = (0:M-1)*Ts; % the time

%Using the given parameteres we reconstruct the input signal
%input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);

prob = zeros(M,1);

if strcmp(runmod.noise_model,'Gauss')
%     disp(sprintf('Calculating Gauss probabilities'));
    if ~isempty(Y_index_others)
        prob(Y_index_others) = .5*(erf( ( Q(Y_index(Y_index_others)+1) - input_signal(Y_index_others))/sigma/sqrt(2)) ...
            - erf(( Q(Y_index(Y_index_others)) - input_signal(Y_index_others))/sigma/sqrt(2)));
    end%if

    if ~isempty(Y_index_0)
        prob(Y_index_0) = 0.5*(1+erf( ( Q(1) - input_signal(Y_index_0))/sigma/sqrt(2)) );
    end%if

    if ~isempty(Y_index_N)
        prob(Y_index_N) = 0.5*(1-erf( ( Q(N) - input_signal(Y_index_N))/sigma/sqrt(2)));
    end%if

elseif strcmp(runmod.noise_model,'Laplace')
%     disp(sprintf('Calculating Laplace probabilities'));
    if ~isempty(Y_index_others)
    Tlm  = Q(Y_index(Y_index_others)+1) - input_signal(Y_index_others);
    Tl1m = Q(Y_index(Y_index_others)  ) - input_signal(Y_index_others);
    
    prob(Y_index_others) = -0.5*sign(Tlm) .* (exp( -sign(Tlm).*Tlm/sigma ) - 1 ) ...
        +0.5*sign(Tl1m) .* (exp( -sign(Tl1m).*Tl1m/sigma ) - 1 );

    end%if

    if ~isempty(Y_index_0)
        Tlm  = Q(1) - input_signal(Y_index_0); 
        prob(Y_index_0) = -.5*sign(Tlm) .* (exp( -sign(Tlm).*Tlm/sigma ) - 1 ) + 0.5 ;
    end%if

    if ~isempty(Y_index_N)
        Tl1m = Q(N) - input_signal(Y_index_N);
        prob(Y_index_N) = 0.5 + 0.5*sign(Tl1m) .* (exp( -sign(Tl1m).*Tl1m/sigma ) - 1 );
    end%if 
end%if