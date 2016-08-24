function [A1,A2,DC] = init_sin_est(meas_data,runmod)

freq = get_sine_freq(meas_data);
Ts = get_Ts(meas_data);
%discrete_freq = freq*Ts;

if nargin < 2
    runmod.verbose = 3;
end %if


% generates initial estimates of the sine wave

[y,tt] = get_measured_data(meas_data);    

M = length(y);
A = ones(M,3);

if isempty(tt)
    tt = 0:M-1;
    tt = tt.';
end%if

tt = tt*Ts;

%tt = (0:M-1).';
%A(:,1) = sin(2*pi*discrete_freq*tt);
%A(:,2) = cos(2*pi*discrete_freq*tt);
A = ones(M,3);
A(:,1) = sin(2*pi*freq*tt);
A(:,2) = cos(2*pi*freq*tt);

P = A\(y);
A1 = P(1);
A2 = P(2);
DC = P(3);


%%%%%%%%%%%%%%%%%%%%%%%% Plot/Debug ASAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;

% A=sqrt(P(1)^2+P(2)^2);
% 
% if      P(2)>0
%     phi=atan(-P(1)/P(2));
% elseif P(2)<0
%     phi=atan(-P(1)/P(2))+pi;
% else %x0(1)==0
%     if P(1)<0, phi=pi/2;
%     else phi=3*pi/2;
%     end
% end
% 
% yy1 = A1*sin(2*pi*freq*tt')+A2*cos(2*pi*freq*tt')+DC;  
% yy2 = A*cos(2*pi*freq*tt'+phi)+DC;
% 
% [X, Rn, Q]=sfit3(y,tt,1/Ts,freq);
% 
% yy_sfit3a = Q.A*cos(2*pi*freq*tt')+Q.B*sin(2*pi*freq*tt')+Q.C;  
% yy_sfit3b = X.A*cos(2*pi*X.f*tt'+X.phi*pi/180)+X.DC;
% 
% Q
% P
% Q.A-P(1)
% Q.B-P(2)
% Q.C-P(3)
% 
% %figure;
% plot(yy_sfit3a,'k'); %hold;
% plot(yy1,'g--');
% %plot(yy2,'r--');
% 
% figure;
% delta1 = yy1-yy_sfit3a;
% delta2 = yy1-yy_sfit3b;
% delta3 = yy2-yy_sfit3a;
% delta4 = yy2-yy_sfit3b;
% delta5 = yy1-yy2;
% delta6 = yy_sfit3a-yy_sfit3b;
% 
% sum(delta1)
% sum(delta4)
% sum(delta2+delta3)
% sum(delta5-delta6)
% 
% plot(delta5); hold;
% plot(-delta6,'r');
% %plot(delta2,'g--');
% %plot(delta3,'c--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot(Y_index);
% hold on;
% plot( A1 *sin(2*pi*discrete_freq*tt)+A2 * cos(2*pi*discrete_freq*tt) + D,'g');

%%%%%%%%%%%%%%%% Verbose 
if runmod.verbose > 1
    disp('Initial sine estimation done.');
end%if

if runmod.verbose > 2
    disp('-----------------------------------------');
end%if

