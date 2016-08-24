function sigma_est = init_sigma_est(Tl,meas_data,A_init,B_init,C_init,runmod)

if nargin < 6
    runmod.verbose = 3;
end %if

y = get_measured_data(meas_data);
input_signal = reconstruct_input_signal(meas_data,A_init,B_init,C_init);
qd = quant(Tl,input_signal)+1;
s = qd-y;
sa = d2c_amp(meas_data,s);

sigma_est = sqrt(var(sa));
sigma_est = sigma_est*1;

%%%%%%%%%%%%%%%% Verbose 
if runmod.verbose > 1
    disp('Initial variance estimation done.');
end%if

if runmod.verbose > 2
    disp('-----------------------------------------');
end%if

