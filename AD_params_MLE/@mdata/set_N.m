function obj = set_N(obj,Nin)


if ~(Nin > 0)
    error('Number of codes must be positive');
end%if

obj.N_bit = log2(Nin);
obj.N = Nin;

