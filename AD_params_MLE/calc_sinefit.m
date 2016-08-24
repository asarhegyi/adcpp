function [X, Rn]=calc_sinefit(data, sample_rate, start_fr)

%CALC_SINEFIT  Four parameter fit of a sine wave to measured data based on sfit4
%
%       [X, Rn]=calc_sinefit(data,sample_rate,start_fr)
%       The function fits a sinewave onto the input data using 
%       the four parameter fitting described in the IEEE 1240-2000
%       standard.
%       It plots the original samples, the fit and the residual vector.
%
%       Input arguments:
%         data: vector of measured samples
%         sample_rate - sampling rate   (numeric)
%         start_fr    - input frequency in Hz (used as a starting point for
%         iteration)
%
%       Output arguments:
%       X: structure of  DC + A * cos(w*tvect + phi)
%         DC   - DC level
%         A    - amplitude of the fitted cosine
%         f    - frequency of the fitted cosine
%         phi  - phase of the fitted cosine
%         erms - RMS error
%       Rn: the residual vector, the difference between the data and the fitted sine-wave


% Written by Zoltán Tamás Bilau, modified by Janos Markus
% $Original Id: sfit4.m,v 3.0 2004/04/19 11:20:09 markus Exp $
% Copyright (c) 2001-2007 by Istvan Kollar and Janos Markus
% All rights reserved.

if ~isnumeric(data), error('Sample vector is not numeric'); end
if length(size(data))>2, error('Dimension of the data vector is more than 2');	end
if min(size(data))>1, error('The data array is not a vector'); end

N=length(data);
if ~all(finite(data)), error('Not all samples are finite'); end
if any(imag(data))~=0, error('Not all samples are real');end
if N<4, error('Less than four samples'); end

%if nargin>1, sample_rate=varargin{1}; end
%if nargin>2, start_fr=varargin{2}; end

mode='dft';
maxCycN=30;
showiter=0;
werror=2*pi*1e-6;


if ~isnumeric(sample_rate),
    error('sample_rate is not numeric');
end
if ~isequal(size(sample_rate),[1,1]),
    error('sample_rate is not a scalar');
end
if sample_rate<=0,
    error('sample_rate is not positive');
end

Ts=1/sample_rate;
time=[0:N-1]*Ts;
unisamp=1; %uniform sampling

data=data(:); %Force column vector

if ~isequal(size(start_fr),[1,1]),
    error('start_fr is not a scalar');
end
if isnan(start_fr) | ~isfinite(start_fr) | start_fr<0
    error('dfmin is not a positive number')
end
w=2*pi*start_fr;


% three parameters (step 0):
D0=[cos(w*time);sin(w*time);ones(1,N)]';
x0=D0\data;
x0(4)=0;

iQ=0; icyc=0;
while 1
    icyc=icyc+1;
    
    D0=[cos(w*time);sin(w*time);ones(1,N);-x0(1)*time.*sin(w*time)+x0(2)*time.*cos(w*time)]';
    
    x0=D0\data;
    if abs(x0(4))<0.05*w, w=w+x0(4); %adjust w
    else w=w+.05*x0(4);
    end
    if abs(x0(4))<werror/Ts break; end
    if icyc>=maxCycN break; end
end


Rn=data-[x0(1)*cos(w*time)+x0(2)*sin(w*time)+x0(3)]';

figure(1), 
subplot(2,1,1)
stairs(time,data,'b');
hold on;
stairs(time,D0(:,1:3)*x0(1:3),'r');
hold off;
legend('Original','Fitted');
xlabel('Time [s]');
ylabel('Amplitude [LSB]')
axis([time(1), time(end), 0, 4100]);
grid on

subplot(2,1,2)
stairs(time,Rn)
legend('Residual');
xlabel('Time [s]');
ylabel('Amplitude [LSB]')
grid on
axis([time(1), time(end), -1.2*max(abs(Rn)), 1.2*max(abs(Rn))]);

%Compose output:
X.DC=x0(3);
X.A=sqrt(x0(1)^2+x0(2)^2);
X.f=w/2/pi;
if      x0(1)>0
    X.phi=atan(-x0(2)/x0(1));
elseif x0(1)<0
    X.phi=atan(-x0(2)/x0(1))+pi;
else %x0(1)==0
    if x0(2)<0, X.phi=pi/2;
    else X.phi=3*pi/2;
    end
end
X.phi=X.phi*(180/pi);
X.erms=sqrt(1/N*Rn'*Rn);

%
%End of calc_sinefit
