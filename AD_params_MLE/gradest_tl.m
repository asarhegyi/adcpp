function [P_out,t_init_tl_out] = gradest_tl(P_in,Y_index,Ts,freq,meas_data,t_init_tl,runmod)

if runmod.verbose >= 1, disp('Backtracking w.r.t transition levels');end;
%default value
if ~isfield(runmod,'tl_alg')
    runmod.tl_alg = 'full';
end%if

% extract varaibles from input data
[A,B,C,tl,sigma] = split_parameter_vector(P_in);
P = P_in;
N = length(tl);
M = length(Y_index); %number of samples
cf_params.Ts = Ts;
cf_params.Y_index = Y_index;
cf_params.freq = freq;
cf_params.meas_data = meas_data;

P_prev_temp = P;

if strcmp(runmod.tl_alg,'full')
     %------------------------------------------------------------------
    
    % Every transition level estimated separately from each others.
    % If the bit number is high (say 8), the calculation could be very
    % slow.
    
    for ii=1:length(tl)
        % report
        if runmod.verbose >= 2
            if length(tl) > 4
                if length(tl) < 10
                    if ii == ceil(length(tl)/4)
                        disp('Progress: 25 %');
                    elseif ii == ceil(length(tl)/2)
                        disp('Progress: 50 %');
                    elseif ii == ceil(3*length(tl)/4)
                        disp('Progress: 75 %');
                    end%if ii == ceil(length(tl)/4)
                else %if length(tl) < 10
                    if ii == ceil(length(tl)/10)
                        disp('Progress: 10 %');
                    elseif ii == ceil(2*length(tl)/10)
                        disp('Progress: 20 %');
                    elseif ii == ceil(3*length(tl)/10)
                        disp('Progress: 30 %');
                    elseif ii == ceil(4*length(tl)/10)
                        disp('Progress: 40 %');
                    elseif ii == ceil(5*length(tl)/10)
                        disp('Progress: 50 %');
                    elseif ii == ceil(6*length(tl)/10)
                        disp('Progress: 60 %');
                    elseif ii == ceil(7*length(tl)/10)
                        disp('Progress: 70 %');
                    elseif ii == ceil(8*length(tl)/10)
                        disp('Progress: 80 %');
                    elseif ii == ceil(9*length(tl)/10)
                        disp('Progress: 90 %');
                    end%if ii == ceil(length(tl)/10)                    
                end%if length(tl) < 10
                    
            end%if length(tl) > 4
        end %if runmod.verbose >= 2, 
        
        tl_grad = zeros(N,1); %D        
        tl_grad(ii) = gradient_wrt_tl(ii,A,B,C,tl,sigma,meas_data,runmod);

        P_grad_tl = tl_part_of_gradient_vector(tl_grad,1:length(tl),length(P)); %D

        if runmod.verbose > 2, disp(sprintf('Backtracking w.r.t tl(%d)',ii));end; %D
        
        P_temp = P;%D
        [P,stat_tl] = backtracking_line_search(runmod,P,P_grad_tl,-P_grad_tl,cf_params,[],[],t_init_tl(ii));%D
        if stat_tl.number_of_inf > 10
            t_init_tl(ii) = t_init_tl(ii)/10;
        elseif stat_tl.number_of_inf > 5
            t_init_tl(ii) = t_init_tl(ii)/2;
        elseif stat_tl.number_of_inf == 0
            t_init_tl(ii) = t_init_tl(ii)*1.2;
        end
        
        if runmod.verbose > 2, disp(sprintf('delta tl(%d) = %f',ii,num2lsb(meas_data,norm(P-P_temp))));end;%D
        [A,B,C,tl,sigma] = split_parameter_vector(P);%D
     
        input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
        prob = calculate_probabilities(tl,sigma,Y_index,input_signal,runmod);% size = M x 1
        if any(~(prob >0))
            disp('**Problem: Invalid state.');
        end%if
        
    end%for ii
    %--------------------------------------------------------------------    
elseif strcmp(runmod.tl_alg,'groupped')
    % The gradient vector contains more than one coordinates
    % The variance guides how much transition levels can be groupped.
    
    step_index = round(5*num2lsb(meas_data,sigma));
    if step_index == 0
        %it means that the variance is extremly small
        step_index = 1;
    end%if
    
    if step_index > N-1
        %extreme big variance, or very small bit number
        step_index = N;
    end%if
    
%     if step_index > 1
%         disp('Hello');
%     end
    for jj=1:step_index
        
        
        if runmod.verbose >= 2
            if step_index > 4
                if step_index < 10
                    if jj == ceil(step_index/4)
                        disp('Progress: 25 %');
                    elseif jj == ceil(step_index/2)
                        disp('Progress: 50 %');
                    elseif jj == ceil(3*step_index/4)
                        disp('Progress: 75 %');
                    end%if ii == ceil(length(step_index)/4)
                else %if length(step_index) < 10
                    if jj == ceil((step_index)/10)
                        disp('Progress: 10 %');
                    elseif jj == ceil(2*(step_index)/10)
                        disp('Progress: 20 %');
                    elseif jj == ceil(3*(step_index)/10)
                        disp('Progress: 30 %');
                    elseif jj == ceil(4*(step_index)/10)
                        disp('Progress: 40 %');
                    elseif jj == ceil(5*(step_index)/10)
                        disp('Progress: 50 %');
                    elseif jj == ceil(6*(step_index)/10)
                        disp('Progress: 60 %');
                    elseif jj == ceil(7*(step_index)/10)
                        disp('Progress: 70 %');
                    elseif jj == ceil(8*(step_index)/10)
                        disp('Progress: 80 %');
                    elseif jj == ceil(9*(step_index)/10)
                        disp('Progress: 90 %');
                    end%if ii == ceil(length(step_index)/10)                    
                end%if length(step_index) < 10
                    
            end%if length(step_index) > 4
        end %if runmod.verbose >= 2, 
        
        current_index_set = jj:step_index:N;
        tl_grad = zeros(N,1);         
        for ii=1:length(current_index_set)
            tl_grad(current_index_set(ii)) = gradient_wrt_tl(current_index_set(ii),A,B,C,tl,sigma,meas_data,runmod);
        end%for ii
        
        P_grad_tl = tl_part_of_gradient_vector(tl_grad,1:length(tl),length(P)); %D
        P_temp = P;%D
        
        t_init_tl_temp = mean(t_init_tl(current_index_set));
        [P,stat_tl] = backtracking_line_search(runmod,P,P_grad_tl,-P_grad_tl,cf_params,[],[],t_init_tl_temp);%D
        
        if stat_tl.number_of_inf > 10
            t_init_tl(current_index_set) = t_init_tl(current_index_set)/10;
        elseif stat_tl.number_of_inf > 5
            t_init_tl(current_index_set) = t_init_tl(current_index_set)/2;
        elseif stat_tl.number_of_inf == 0
            t_init_tl(current_index_set) = t_init_tl(current_index_set)*1.2;
        end
        
        if runmod.verbose > 2, disp(sprintf('delta tl(%d) = %f',ii,num2lsb(meas_data,norm(P-P_temp))));end;%D
        [A,B,C,tl,sigma] = split_parameter_vector(P);%D
     
        %input_signal = reconstruct_input_signal(A,B,C,M,Ts,freq);
        input_signal = reconstruct_input_signal(meas_data,A,B,C);
        prob = calculate_probabilities(tl,sigma,Y_index,input_signal,runmod);% size = M x 1
        if any(~(prob >0))
            disp('**Problem: Invalid state.');
        end%if
        
    end%foir jj
    %--------------------------------------------------------------------
else
    error('Invalid algorithm (runmod.tl_alg)');
end%strcmp(runmod.tl_alg,'full')


if runmod.verbose >= 1, disp(sprintf('delta tl = %f',num2lsb(meas_data,norm(P-P_prev_temp))));end;%D

% calculate outputs
P_out = P;
t_init_tl_out = t_init_tl;
