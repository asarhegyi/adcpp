function obj = set_Nbit(obj,Nb)

if Nb ~= floor(Nb)
    error('Bit number must be an integer value');
end%if

if ~(Nb > 0)
    error('Bit number must be positive');
end%if

obj.N_bit = Nb;
obj.N = 2^Nb;

