function [X, Rn, Q]=sfit4(data, varargin)
%SFIT4  IEEE-STD-1241 Standard four parameter fit of a sine wave to measured data
%
%       [X, Rn, Q]=sfit4(data,time,sample_rate,start_fr)
%
%       Input arguments:
%         data: vector of measured samples
%       
%       Optional parameters (some optional parameters can be either numeric or character string):
%         time: vector of time spacing (numeric)
%         sample_rate: sampling rate   (numeric)
%         start_fr - starting frequency in Hz (numeric or char string) (Default is determined by DFT or IpFFT)
%
%       Output arguments:
%       X: structure of  DC + A * cos(w*tvect + phi)
%         DC   - DC level
%         A    - amplitude of the fitted cosine
%         f    - frequency of the fitted cosine
%         phi  - phase of the fitted cosine
%         erms - RMS error
%       Rn: the residual vector, the difference between the data and the fitted sine-wave
%       Q: structure array of iteration steps
%
% See also SFIT3, SFIT4IMP


% Written by Zoltan Tamas Bilau, modified by Janos Markus
% $Id: sfit4.m,v 3.0 2004/04/19 11:20:09 markus Exp $
% Copyright (c) 2001-2004 by Istvan Kollar and Janos Markus
% All rights reserved.

if ~isnumeric(data), error('Sample vector is not numeric'); end
if length(size(data))>2, error('Dimension of the data vector is more than 2');    end
if min(size(data))>1, error('The data array is not a vector'); end

N=length(data);
if ~all(isfinite(data)), error('Not all samples are finite'); end
if any(imag(data))~=0, error('Not all samples are real');end
if N<4, error('Less than four samples'); end
mode='';

if nargin>1, time=varargin{1}; end
if nargin>2, sample_rate=varargin{2}; end
if nargin>3, start_fr=varargin{3}; end

mode='dft';
maxCycN=30;
showiter=0;
werror=2*pi*1e-6;


unisamp=0; %uniform sampling
if exist('sample_rate','var') 
    if ~isnumeric(sample_rate), 
        warning('sample_rate is not numeric'); 
    end
    if ~isequal(size(sample_rate),[1,1]), 
        warning('sample_rate is not a scalar'); 
    end
    if sample_rate<=0, 
        warning('sample_rate is not positive'); 
    end
end


if exist('time','var')
    if isempty(time)
        if ~exist('Ts'), Ts=1; sample_rate=1; end
        time=[0:N-1]*Ts; % Ts default
        unisamp=1;
    end
    if ~isnumeric(time), error('time is not numeric'); end
    if ~isequal(size(data),size(time)), error('sizes of data and time differ'), end
    time=time(:)';
    %Unisamp is 0 even if time is given and has almost equal distances
    unisamp=all(diff(diff(time))<max(diff(time))*100*eps);
    if ~exist('Ts','var')
        Ts=median(diff(time));
        sample_rate=1/Ts;
    end
else
    if ~exist('Ts'), Ts=1; sample_rate=1; end
    time=[0:N-1]*Ts; % Ts default
    unisamp=1;
end

data=data(:); %Force column vector

%If startW is given
w=0;
if exist('start_fr','var')
    %evaluate parameter
    if isstr(start_fr)
        if isempty(start_fr) | all(isspace(start_fr))
            start_fr=0;
        else
            try
                start_fr=eval(start_fr);
            catch
                warning('Cannot evaluate start_fr, initial frequency is determined by DFT')
                start_fr=0;
            end
        end
    elseif ~isnumeric(start_fr) 
        warning('Cannot interpret start_fr, initial frequency is determined by DFT')
        start_fr=0;        
    end
    
    if ~isequal(size(start_fr),[1,1]), 
        warning('start_fr is not a scalar, initial frequency is determined by DFT'); 
        start_fr=0;
    end
    
    if isnan(start_fr) | ~isfinite(start_fr) | start_fr<0
        warning('dfmin<0, initial frequency is determined by DFT')
        start_fr=0;
    end
    w=2*pi*start_fr;
end

if ~w %w is 0, no initial frequency is given
    %Initial radian frequency: look at maximum of PSD
    
    if ~unisamp
        Nt=round((max(time)-min(time)+Ts)/Ts);
        tvecti=round(time/Ts);
        
        %interpolation should be preferred if there are missing data
        tvecti=tvecti-min(tvecti)+1; 
        
        ym=zeros(Nt,1);
        ym(tvecti)=data;
        ip=[1:Nt]'; ip(tvecti)=[];
        for ii=1:length(ip)
          if ym(min(end,ip(ii)+1))~=0, ym(ip(ii))=ym(min(end,ip(ii)+1));
          elseif ym(max(1,ip(ii)-1))~=0, ym(ip(ii))=ym(max(1,ip(ii)-1));
          end
        end
        for ii=1:length(ip):-1:1
          if ym(ip(ii))==0
            if ym(min(end,ip(ii)+1))~=0, ym(ip(ii))=ym(min(end,ip(ii)+1));
            elseif ym(max(1,ip(ii)-1))~=0, ym(ip(ii))=ym(max(1,ip(ii)-1));
            end
          end
        end
        figure(1), plot(time,data,'.','markersize',0.4)
        ylim=get(gca,'ylim');
        ylim(2)=pow2(ceil(log2(ylim(2))));
        if ylim(1)>0
          ylim(1)=pow2(floor(log2(ylim(1))));
          if ylim(1)<8, ylim(1)=0; end
        elseif ylim(1)<0
          ylim(1)=-pow2(ceil(log2(abs(ylim(1)))));
        end
        set(gca,'ylim',ylim)
        title('Processed samples'), zoom on
        F=abs(fft(ym));
        [Mfft,w]=max(F(2:round(Nt/2)));
        %figure(99), plot(ym), xlabel(sprintf('Index of max: %.0f',w)), shg
        w=2*pi*w/(Nt*Ts);
    else
        Fc=fft(data);
        F=abs(Fc);
        [Mfft,w]=max(F(2:round(N/2)));
        
        if strcmpi(mode,'ipfft')
            if w>1
                %calculating the 2 points, between them the estimated frequency is
                if F(w-1)>F(w+1) w=w-1; end
            end
            
            n=2*pi/N;
            U=real(Fc(w+1));    V=imag(Fc(w+1));
            U1=real(Fc(w+2));  V1=imag(Fc(w+2));
            Kopt=(sin(n*w)*(V1-V)+cos(n*w)*(U1-U))/(U1-U);
            Z1=V*(Kopt-cos(n*w))/sin(n*w)+U;
            Z2=V1*(Kopt-cos(n*(w+1)))/sin(n*(w+1))+U1;
            
            lambda=acos((Z2*cos(n*(w+1))-Z1*cos(n*w))/(Z2-Z1))/n;
            w=lambda;
        end
        w=2*pi*w/(N*Ts);
    end
end


% three parameters (step 0):
D0=[cos(w*time);sin(w*time);ones(1,N)]';
x0=D0\data;
figure(1), hold on, hp=plot(time,D0*x0,'.g','markersize',.4); hold off, shg
x0(4)=0;

iQ=0; icyc=0;
while 1
    icyc=icyc+1;
    if nargout>2 
        iQ=iQ+1;
        Q(iQ).A=x0(1);
        Q(iQ).B=x0(2);
        Q(iQ).C=x0(3);
        Q(iQ).w=w;
    end
    
    D0=[cos(w*time);sin(w*time);ones(1,N);-x0(1)*time.*sin(w*time)+x0(2)*time.*cos(w*time)]';
    
    x0=D0\data;
    if abs(x0(4))<0.05*w, w=w+x0(4); %adjust w
    else w=w+.05*x0(4);
    end
    if showiter
        if icyc==1
            disp(' ');
            disp('   Cyc#   frequency in Hz   Change of frequency              A           B');
            disp('   ----   ---------------   ----------------------------------------------');
        end
        fprintf('   %02d     %1.12g  %+1.25g             %+1.10g  %+1.10g\n',icyc,w/2/pi,Ts*x0(4)/2/pi,x0(1),x0(2));
    end
    
    if abs(x0(4))<werror/Ts break; end
    if icyc>=maxCycN break; end
end

figure(1), delete(hp), hp=plot(time,D0(:,1:3)*x0(1:3),'.r','markersize',.4); hold off, shg
close(1)

if nargout<2
    iQ=iQ+1;
    Q(iQ).A=x0(1);
    Q(iQ).B=x0(2);
    Q(iQ).C=x0(3);
    Q(iQ).w=w;
end

Rn=data-[x0(1)*cos(w*time)+x0(2)*sin(w*time)+x0(3)]';

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
%End of sfit4
