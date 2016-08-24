function [A,B,C,Tl,sigma] = ML_est_ABCN(meas_data,runmod)

% This function calculates the maximum likelihood estimation.
%
% Input arguments:
%       meas_data : type of class mdata. Describes the measurement. For
%           example tha maximal analog transition level etc.
%       runmod : run modifier
%
% Output arguments:
%       A,B,C : result of the estimation w.r.t. parameters of the sine
%           wave. The assumed signal is
%           y(t) = A*sin(2*pi*freq*t)+B*cos(2*pi*freq*t)+C
%       Tl : the transition levels
%       sigma : deviation of the additive Gaussian noise

if nargin < 2
    %in this case we have no runmod
    runmod.verbose = 3;
    runmod.eps_limit = 0.01;
    runmod.noise_model = 'Gauss';
end%if;

if runmod.verbose > 1
    disp('-----------------------------------------');
    disp('Estimation algorithm started.');
end%if

freq = get_sine_freq(meas_data);
if isempty(freq)
    error('Four parameter fitting has not been implemented yet');
end%if


%Inital estimate of the transition levels (Tl).
%These are the ideal transition levels. It is important to have an initial estimate.
Tl_ideal = get_Tl_of_an_ideal_quantizer(meas_data);

%Y_index = Y_index-min(Y_index);
%COMMENTED THIS OUT because I think this is also done in remove_missing_code.m.  Attila 
%meas_data = corrigate_minium_value_of_samples(meas_data); %might not need this compare with remove_missing_code 

%It is very important to remove missing code.
%To do this the indexes are modified and the transition levels are also
%corrected. Tl will replace Tl_ideal. It can be done because after this
%every parameter is in the analog domain.
Y_index = get_measured_data(meas_data);     %get the indexes
Y_index_save = Y_index;
[Y_index,Tl] = remove_missing_code(Y_index_save,Tl_ideal);
meas_data = update_meas_data(meas_data,Y_index);

%Initial estimation of the of the sine wave parameters. 
%It is done in the digital domain.
%During sine fitting the transition levels are not used.
[Ad,Bd,Cd] = init_sin_est(meas_data,runmod);

A = d2c_amp(meas_data,Ad);
B = d2c_amp(meas_data,Bd);
C = d2c(meas_data,Cd);

sigma = init_sigma_est(Tl,meas_data,A,B,C,runmod);

%%%%%%%%%%%%%%%%%%%%%%%% Plot/Debug ASAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if runmod.verbose > 2

    %Plot initial estimate of the input sine wave
    M = length(Y_index);                        %number of samples
    plot(0:M-1,Y_index,'g');

    freq = get_sine_freq(meas_data);
    Ts = get_Ts(meas_data);
    [y,tt] = get_measured_data(meas_data);
    tt = tt*Ts;

    yy1 = Ad*sin(2*pi*freq*tt')+Bd*cos(2*pi*freq*tt')+Cd;
    yy2 = A*sin(2*pi*freq*tt')+B*cos(2*pi*freq*tt')+C;

    plot(yy1,'k--');
    plot(20*yy2,'b--');
    
end%if
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(runmod.sine,'AphiDC')
%Convert A,B, and C into A, phi, and DC

    if B>0
        phi=atan(-A/B);
    elseif B<0
        phi=atan(-A/B)+pi;
    else %A==0
        if A<0, phi=pi/2;
        else phi=3*pi/2;
        end
    end
    A=sqrt(A^2+B^2);
    DC = C;
    x0=[A,phi,DC];
elseif strcmp(runmod.sine,'ABC')
    x0=[A,B,C];
end%if


%%%%%%%%%%%%%%%%%%%%%%%%% This is the optimization part %%%%%%%%%%%%%%%%%%%%%%%%% 

if strcmp(runmod.search_alg,'fminsearch')
    
    if runmod.verbose > 2
        options=optimset('Display','iter','FunValCheck','on');
    else
        options=optimset('Display','final','FunValCheck','on');
    end%if
    
    if runmod.noise == 1
        x0=[x0,sigma];
        [xopt,fopt,exitflag]=fminsearch(@(x) Likelihood_wrt_ABCN(x,Tl,meas_data,runmod),x0,options);
    else
        [xopt,fopt,exitflag]=fminsearch(@(x) Likelihood_wrt_ABC(x,Tl,sigma,meas_data,runmod),x0,options);
    end%if

elseif strcmp(runmod.search_alg,'fminunc')

    if runmod.verbose > 2
        oldopts=optimset('fminunc');
        options=optimset(oldopts,'LargeScale','off','GradObj','on','HessUpdate','bfgs','display','iter');
        % options=optimset(oldopts,'LargeScale','off','GradObj','on','HessUpdate','bfgs','display','none','showstatus','iterplus');
    else
        oldopts=optimset('fminunc');
        options=optimset(oldopts,'LargeScale','on','GradObj','on','HessUpdate','bfgs');
    end%if
    
    if runmod.noise == 1
        x0=[x0,sigma];
        [xopt,fopt,exitflag]=fminunc(@(x) Likelihood_wrt_ABCN(x,Tl,meas_data,runmod),x0,options);
    else
        [xopt,fopt,exitflag]=fminunc(@(x) Likelihood_wrt_ABC(x,Tl,sigma,meas_data,runmod),x0,options);
    end%if
    
end%if

A = xopt(1);
B = xopt(2);
C = xopt(3);


if runmod.noise == 1
    sigma = xopt(4);
end%if



% scale1=0.98:0.001:1.02;
% 
% for k=1:length(scale1)
% %    A=scale1(k)*A;
% %    B=scale1(k)*B;
% %    C=scale1(k)*C;
% %    x=[A,B,C];
% %    [L(k),grad] = Likelihood_wrt_ABC(x,Tl,sigma,meas_data)
% 
% %    A=scale1(k)*A;
% %    phi=scale1(k)*phi;
%    DC=scale1(k)*DC;
%    x=[A,phi,DC];
%    [L(k),grad] = Likelihood_wrt_ABC(x,Tl,sigma,meas_data)
% end
% 
%     
%     
% %     x=[A,B,C,sigma];
% %     L(k) = Likelihood_wrt_ABCN(x,Tl,meas_data);
% 
% 
% plot(L);
