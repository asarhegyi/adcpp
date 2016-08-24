function Tl = get_Tl_of_an_ideal_quantizer(obj)

%N = 2.^(obj.N_bit);
N = obj.N;
Q = (obj.V_max-obj.V_min)/N;

Tl = obj.V_min-Q/2 + Q*[1:N-1];
Tl = Tl.';
