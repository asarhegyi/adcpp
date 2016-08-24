function [Y_index_mod,Tl] = remove_missing_code(Y_index_in,Tl_ideal)

if min(Y_index_in) ~= 0
    Y_index_in = Y_index_in-min(Y_index_in);
end%if

code_distribution = nof_symbol(Y_index_in,max(Y_index_in));

N_ideal = length(Tl_ideal)+1;

if isempty(find(code_distribution == 0))
    % no missing code
    Tl = Tl_ideal;
    Y_index_mod  = Y_index_in;
else
    % we have at least one missing code
    Tl = Tl_ideal(find(code_distribution(2:end) ~= 0)); %remove those which are not in the samples
    %implicitly we know that code_distribution[1] ~= 0 and
    %code_distribution[max_index] ~= 0
    
    S = sum(code_distribution ~= 0);% nof symbols
    
    % prepare remap
    index_map = zeros(N_ideal,1);
    index_map( find(code_distribution ~= 0)) = 1:(S);
    
    Y_index_mod = index_map(Y_index_in+1)-1 ;
end%if




