%% MATLAB
function out = plot_stimulus(args)


% This script prints figures to visualize the stimulus generated
% by gen_stimulus.cpp for verification and documentation purposes.
%
%        out = plot_stimulus(args)
%
%       Input arguments:
%         args.arg1: number of quantizer bits
%
%       This program is public domain. It is available through
%       https://github.com/asarhegyi/

% $Id:  $
% Copyright (c) 2015-2016 by Attila Sarhegyi
% All rights reserved.


%close all;

qbits=args.arg1;


if qbits == 3
    compFactor = 0.90;
elseif (qbits == 4)
    compFactor = 0.97;
else
    compFactor = 1;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% plot noisy signal path with ideal quantizer %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data to verify signal path
data = load('data.dat');
noisyData = load('noisyData.dat');
quantizedData = load('quantizedData.dat');

h0 = figure;
t = 1:length(data);
p0 = plot(t,noisyData,'r.');
hold;
[ax,p1,p2] = plotyy(t,data,t,quantizedData,'plot','stairs');
xlabel(ax(1),'Samples') % label x-axis
ylabel(ax(1),'Analog Input [V]') % label left y-axis
ylabel(ax(2),'Quantized Output') % label right y-axis

%set(ax(1),'xlim',[0 600]);
y1max = 2^qbits/(2^qbits-1)*max(data)*compFactor;
y2max = 2^qbits;
set(ax(1),'ylim',[min(data) y1max]);
set(ax(2),'ylim',[min(quantizedData) y2max]);
set(ax(1),'ytick',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
if qbits == 3
    set(ax(2),'ytick',[0 2 4 7]);
elseif (qbits == 4)
    set(ax(2),'ytick',[0 5 10 15]);
end

set(p1,'LineStyle','--');
%set(p0,'LineWidth',2);
set(p1,'LineWidth',1);
set(p2,'LineWidth',1);

title('Noisy Signal Path');
legend('Noisy Stimulus','Noise-free Stimulus','Quantized Output','Location','Best');
grid(ax(1),'on');
qfigset;
%print -deps noisy_signal_path.eps;
print -depsc2 noisy_signal_path.eps;
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plot noisy CTL path with noise free singal %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data to verify CTL path
idealTl = load('quantizedData_idealTl.dat');
noisyTl = load('quantizedData_noisyTl.dat');

h1 = figure;
t = 1:length(data);
p0 = stairs(noisyTl,'r');
hold;
[ax,p1,p2] = plotyy(t,idealTl,t,data,'stairs','plot');
xlabel(ax(1),'Samples') % label x-axis
ylabel(ax(1),'Quantized Output') % label left y-axis
ylabel(ax(2),'Analog Input [V]') % label right y-axis

%set(ax(1),'xlim',[0 600]);
%set(ax(2),'xlim',[0 600]);
y1max = 2^qbits;
y2max = 2^qbits/(2^qbits-1)*max(data)*compFactor;
set(ax(1),'ylim',[min(idealTl) y1max]);
set(ax(2),'ylim',[min(data) y2max]);
if qbits == 3
    set(ax(1),'ytick',[0 2 4 7]);
elseif (qbits == 4)
    set(ax(1),'ytick',[0 5 10 15]);
end    
set(ax(2),'ytick',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);

set(p1,'LineStyle','--');
set(p2,'LineStyle','-.');
set(p0,'LineWidth',2);
set(p1,'LineWidth',1);
set(p2,'LineWidth',1);

title('Noisy CTL Path');
legend('Noisy Quantizer','Ideal Quantizer','Noise-free stimulus','Location','Best');
grid(ax(2),'on');
qfigset;
%print -deps noisy_ctl_path.eps;
print -depsc2 noisy_ctl_path.eps;
hold off;








