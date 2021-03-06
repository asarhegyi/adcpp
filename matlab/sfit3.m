
function [X, Rn, Q]=sfit3(data, varargin)


%SFIT3  IEEE-STD-1241 standard three parameter fit of a sine wave to measured data
%
%       [X, Rn, Q]=sfit4(data,time,sample_rate,input_fr)
%
%       Input arguments:
%         data: vector of measured samples
%         time: vector of time spacing (numeric)
%         sample_rate: sampling rate   (numeric)
%         input_fr - starting frequency in Hz (numeric or char string)
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
% See also SFIT4, SFIT4IMP


% Written by Zolt�n Tam�s Bilau, modified by Janos Markus
% $Id: sfit3.m,v 3.0 2004/04/19 11:20:09 markus Exp $
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
if nargin>3, input_fr=varargin{3}; end


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
        sample_rate=1;
    end
    Ts=1/sample_rate;
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
    
    %Unisamp is 0 even if time is given and generated by
    %time=1/sample_rate*(0:length(data)-1);
    
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

%input_fr must be given
w=0;
if exist('input_fr','var')
    %evaluate parameter
    if isstr(input_fr)
        if isempty(input_fr) | all(isspace(input_fr))
            warning('input_fr is empty, default (1000) is used')
            input_fr=1000;
        else
            try
                input_fr=eval(input_fr);
            catch
                warning('Cannot evaluate input_fr, default (1000) is used')
                input_fr=1000;
            end
        end
    elseif ~isnumeric(input_fr) 
        warning('Cannot interpret input_fr, default (1000) is used')
        input_fr=1000;        
    end
    
    if ~isequal(size(input_fr),[1,1]), 
        warning('input_fr is not a scalar, default (1000) is used'); 
        input_fr=1000;
    end
    
    if isnan(input_fr) | ~isfinite(input_fr) | input_fr<0
        warning('dfmin<0, default (1000) is used')
        input_fr=1000;
    end
    w=2*pi*input_fr;
else
    warning('cannot find input_fr, default (1000) is used')
    input_fr=1000;
    w=2*pi*input_fr;
end


% three parameter LS fit:
D0=[cos(w*time);sin(w*time);ones(1,N)]';
x0=D0\data;
x0(4)=0;

Rn=data-[x0(1)*cos(w*time)+x0(2)*sin(w*time)+x0(3)]';

%Compose output:
X.DC=x0(3);
X.A=sqrt(x0(1)^2+x0(2)^2);
X.f=input_fr;
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
%End of sfit3
