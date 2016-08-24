%cd 'C:\Documents and Settings\asar\My Documents\PhD_thesis\MATLAB\'

close all;
clear all;
clc;

Vmax =  2;
Vmin = -2;
N = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%% Input signal from the Monte Carlo sim

load MC_results_Gauss.txt
Gauss_Sine  = MC_results_Gauss(:,1:3);
Gauss_Tk    = MC_results_Gauss(:,4:2^N+2);
Gauss_Noise = MC_results_Gauss(:,2^N+3);
Gauss_Time  = MC_results_Gauss(:,2^N+4);

% load MC_results_Laplace.txt
% Laplace_Sine  = MC_results_Laplace(:,1:3);
% Laplace_Tk    = MC_results_Laplace(:,4:2^N+2);
% Laplace_Noise = MC_results_Laplace(:,2^N+3);
% Laplace_Time  = MC_results_Laplace(:,2^N+4);

% %%%%%%%%%%%%%%%%%%%%%%%%%% Generate the ideal code trnasition levels
% 
% M = 2^N;
% Q = (Vmax-Vmin)/M;
% Tk_ideal = Vmin-Q/2 + Q*(1:M-1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Mean of Tk for 1000 Monte Carlo runs
%  
% Gauss_Tk_mean    = mean(Gauss_Tk,1);            %average of the MC runs for the Gauss noise model
% Gauss_delta_Tk   = (Tk_ideal - Gauss_Tk_mean)/Q;
% Laplace_Tk_mean  = mean(Laplace_Tk,1);          %average of the MC runs for the Laplace noise model
% Laplace_delta_Tk = (Tk_ideal - Laplace_Tk_mean)/Q;
% 
% h1 = figure;
% plot(Gauss_delta_Tk,'LineWidth',[2]); hold;
% plot(Laplace_delta_Tk,'.r','LineWidth',[2]);
% grid on;
% hc1 = get(h1,'children');
% set(hc1,'fontsize',16);
% legend('Gaussian noise model','Laplace noise model','Location','Best');
% Title('MC simulation results for Tk');
% xlabel('Transition levels');
% ylabel('Tid-Test [LSB]');
% set(hc1,'xlim',[1 M]);
% qfigset;
% %print -depsc2 MC_mean_Tk_GvsL.eps;
% print -deps MC_mean_Tk_GvsL.eps;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Standard Deviation of Tk for 1000 Monte Carlo runs
% 
% Gauss_Tk_sigma   = std(Gauss_Tk)/Q;
% Laplace_Tk_sigma = std(Laplace_Tk)/Q;
% 
% h2 = figure;
% plot(Gauss_Tk_sigma,'LineWidth',[2]); hold;
% plot(Laplace_Tk_sigma,'.r','LineWidth',[2]);
% grid on;
% hc2 = get(h2,'children');
% set(hc2,'fontsize',16);
% legend('Gaussian noise model','Laplace noise model','Location','Best');
% Title('MC simulation results for Tk');
% xlabel('Transition levels');
% ylabel('std(Test) [LSB]');
% set(hc2,'xlim',[1 M]);
% qfigset;
% %print -depsc2 MC_std_Tk_GvsL.eps;
% print -deps MC_std_Tk_GvsL.eps;


%%%%%%%%%%%%%%%%%%%%%%%%%% Sine Wave 


freq = 1;
dc_level = 0;
amplitude = 2.021;
Ts = 1/10e3;
N_samples = 8*8000;
sigma = .018;

% Input sine wave (noiseless)
tt = 0:N_samples-1;
tt = tt*Ts;
zz = amplitude*sin(2*pi*freq*tt+pi*1.2)+dc_level;

%Compose estimated sine wave for the Gauss noise model 
G0 = Gauss_Sine;
G.DC = G0(:,3);
G.A = sqrt( G0(:,1).^2 + G0(:,2).^2 );
G.f = freq;
for k = 1:length(G0) 
    if      G0(k,2)>0
        G.phi(k)=atan(-G0(k,1)/G0(k,2));
    elseif G0(k,2)<0
        G.phi(k)=atan(-G0(k,1)/G0(k,2))+pi;
    else %G0(1)==0
        if G0(k,1)<0, G.phi(k)=pi/2;
        else G.phi(k)=3*pi/2;
        end
    end
end

Gauss_Delta.A = (amplitude - G.A)/amplitude;
Gauss_Delta.phi = (pi*1.2 - G.phi)/(pi*1.2);
Gauss_Delta.DC = (dc_level - G.DC);

% %Compose estimated sine wave for the Laplace noise model 
% L0 = Laplace_Sine;
% L.DC = L0(:,3);
% L.A = sqrt( L0(:,1).^2 + L0(:,2).^2 );
% L.f = freq;
% for k = 1:length(L0) 
%     if      L0(k,1)>0
%         L.phi(k)=atan(-L0(k,2)/L0(k,1));
%     elseif L0(k,1)<0
%         L.phi(k)=atan(-L0(k,2)/L0(k,1))+pi;
%     else %L0(1)==0
%         if L0(k,2)<0, L.phi(k)=pi/2;
%         else L.phi(k)=3*pi/2;
%         end
%     end
% end
% 
% Laplace_Delta.A = (amplitude - L.A)/amplitude;
% Laplace_Delta.phi = (pi*1.2 - L.phi)/(pi*1.2);
% Laplace_Delta.DC = (dc_level - L.DC);

% This is how you can get back the cosine; Do this later matrix mismatch
% zz_est_Gauss = G.A*cos(2*pi*G.f*tt+G.phi)+G.DC;

% plot estimation error of the amplitude for Gauss
h3 = figure;
hist(Gauss_Delta.A,length(Gauss_Delta.A)/5);
hp3 = findobj(gca,'Type','patch');
set(hp3,'FaceColor','black','EdgeColor','none');
grid on;
hc3 = get(h3,'children');
set(hc3,'fontsize',16);
%legend('Gaussian noise model','Location','Best');
Title('Estimation error of the amplitude');
xlabel('(A-Aest)/A');
ylabel('Frequency');
qfigset;
%print -depsc2 MC_Amp_GvsL.eps;
print -deps MC_Amp_G.eps;


% % plot estimation error of the amplitude for Laplace
% h3 = figure;
% hist(Laplace_Delta.A,length(Laplace_Delta.A)/5);
% hp3 = findobj(gca,'Type','patch');
% set(hp3,'FaceColor','black','EdgeColor','none');
% grid on;
% hc3 = get(h3,'children');
% set(hc3,'fontsize',16);
% %legend('Laplace noise model','Location','Best');
% Title('Estimation error of the amplitude');
% xlabel('(A-Aest)/A');
% ylabel('Frequency');
% qfigset;
% %print -depsc2 MC_Amp_GvsL.eps;
% print -deps MC_Amp_L.eps;

% figure;
% plot(Gauss_Delta.phi); hold; plot(Laplace_Delta.phi,'r');
% figure;
% plot(Gauss_Delta.DC); hold; plot(Laplace_Delta.DC,'r');

% calculcate mean and std of the estimation error of A, phi, and DC 
Gauss_Delta_mean.A = mean((amplitude - G.A)/amplitude);
Gauss_Delta_mean.phi = mean((pi*1.2 - G.phi)/(pi*1.2));
Gauss_Delta_mean.DC = mean((dc_level - G.DC));

Gauss_Delta_std.A = std((amplitude - G.A)/amplitude);
Gauss_Delta_std.phi = std((pi*1.2 - G.phi)/(pi*1.2));
Gauss_Delta_std.DC = std((dc_level - G.DC));

% Laplace_Delta_mean.A = mean((amplitude - L.A)/amplitude);
% Laplace_Delta_mean.phi = mean((pi*1.2 - L.phi)/(pi*1.2));
% Laplace_Delta_mean.DC = mean((dc_level - L.DC));
% 
% Laplace_Delta_std.A = std((amplitude - L.A)/amplitude);
% Laplace_Delta_std.phi = std((pi*1.2 - L.phi)/(pi*1.2));
% Laplace_Delta_std.DC = std((dc_level - L.DC));

disp('##############');
disp(sprintf('The average amplitude estimation error for the Gauss noise model is %.2e with %.2e std',Gauss_Delta_mean.A, Gauss_Delta_std.A));
%disp(sprintf('The average amplitude estimation error for the Laplace noise model is %.2e with %.2e sec std',Laplace_Delta_mean.A, Laplace_Delta_std.A));
disp('##############');
disp(sprintf('The average phase estimation error for the Gauss noise model is %.2e with %.2e std',Gauss_Delta_mean.phi, Gauss_Delta_std.phi));
%disp(sprintf('The average phase estimation error for the Laplace noise model is %.2e with %.2e sec std',Laplace_Delta_mean.phi, Laplace_Delta_std.phi));
disp('##############');
disp(sprintf('The average DC estimation error for the Gauss noise model is %.2e with %.2e std',Gauss_Delta_mean.DC, Gauss_Delta_std.DC));
%disp(sprintf('The average DC estimation error for the Laplace noise model is %.2e with %.2e sec std',Laplace_Delta_mean.DC, Laplace_Delta_std.DC));
disp('##############');


% %%%%%%%%%%%%%%%%%%%%%%%%%% Excecution Time 
% 
% 
% G_time_avg = mean(Gauss_Time);
% G_time_std = std(Gauss_Time');
% 
% L_time_avg = mean(Laplace_Time);
% L_time_std = std(Laplace_Time');
% 
% disp(sprintf('The average excecution time for the Gauss noise model is %.2f sec with %.2f sec std',G_time_avg, G_time_std));
% disp(sprintf('The average excecution time for the Laplace noise model is %.2f sec with %.2f sec std',L_time_avg, L_time_std));
% disp('##############');
% 
% figure;
% hist(Gauss_Time(:,1),length(Gauss_Time)/5);
% title('Histogram of the excecution time for Gauss');
% qfigset;
% 
% figure;
% hist(Laplace_Time(:,1),length(Laplace_Time)/5);
% title('Histogram of the excecution time for Laplace');
% qfigset;