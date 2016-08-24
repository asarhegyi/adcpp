function out_lsb = num2lsb(obj,in_num)

%N = 2.^(obj.N_bit);
N = obj.N;
lsb = (obj.V_max-obj.V_min)/N;

out_lsb = in_num/lsb;
