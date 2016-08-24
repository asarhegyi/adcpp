function code_density_INL_DNL(meas_data,runmod)

if nargin < 2
    %in this case we have no runmod
    runmod.verbose = 3;
    runmod.eps_limit = 0.01;
    runmod.noise_model = 'Gauss';
end%if;

if runmod.verbose > 1
    disp('-----------------------------------------');
    disp('Code Density algorithm started.');
end%if


% filename=input('Enter File Name: ');
% if isempty(filename)
%    filename = 'listing';
% end
% fid=fopen(filename,'r');
% numpt=input('Enter Number of Data Points:  ');
% numbit=input ('Enter ADC Resolution:  ');
% mid_code=input('Enter Mid-Code (Mean Code):  ');


numpt=get_measured_data_length(meas_data);
numbit=8;

[code,time] = get_measured_data(meas_data);     %get the indexes
code_count=zeros(1,2^numbit);

for i=1:size(code),
   code_count(code(i)+1)=code_count(code(i)+1) + 1;
end

if code_count(1) == 0 | code_count(2^numbit) == 0 | ...
  (code_count(1) < code_count(2)) | (code_count(2^numbit-1) > code_count(2^numbit))
   disp('ADC not clipping ... Increase sinewave amplitude!');
   return;
end

A=max((2^numbit)/2,2^numbit-1-(2^numbit)/2); 
vin=(0:2^numbit-1)-(2^numbit)/2;	
sin2ramp=1./(pi*sqrt(A^2*ones(size(vin))-vin.*vin));


while sum(code_count(3:2^numbit-2)) < numpt*sum(sin2ramp(3:2^numbit-2))
  A=A+0.01;
  sin2ramp=1./(pi*sqrt(A^2*ones(size(vin))-vin.*vin));
end


disp('You Have Applied a Sine Wave of (dBFS): '); 
Amplitude=A/(2^numbit/2)
figure;
plot([0:2^numbit-1],code_count,[0:2^numbit-1],sin2ramp*numpt);
title('CODE HISTOGRAM - SINE WAVE');
xlabel('DIGITAL OUTPUT CODE');
ylabel('COUNTS');
axis([0 2^numbit-1 0 max(code_count(2),code_count(2^numbit-1))]);
code_countn=code_count(2:2^numbit-1)./(numpt*sin2ramp(2:2^numbit-1)); 
figure;
plot([1:2^numbit-2],code_countn);
title('CODE HISTOGRAM - NORMALIZED')
xlabel('DIGITAL OUTPUT CODE');
ylabel('NORMALIZED COUNTS');

dnl=code_countn-1;
inl=zeros(size(dnl));
for j=1:size(inl')
   inl(j)=sum(dnl(1:j));
end

%End-Point fit INL
%[p,S]=polyfit([1,2^numbit-2],[inl(1),inl(2^numbit-2)],1);

%Best-straight-line fit INL
[p,S]=polyfit([1:2^numbit-2],inl,1);
inl=inl-p(1)*[1:2^numbit-2]-p(2);

disp('End Points Eliminated for DNL and INL Calculations');
figure;
plot([1:2^numbit-2],dnl);
grid on;
title('DIFFERENTIAL NONLINEARITY vs. DIGITAL OUTPUT CODE');
xlabel('DIGITAL OUTPUT CODE');
ylabel('DNL (LSB)');
figure;
plot([1:2^numbit-2],inl);
grid on;
title('INTEGRAL NONLINEARITY vs. DIGITAL OUTPUT CODE');
xlabel('DIGITAL OUTPUT CODE');
ylabel('INL(LSB)');  

disp('Done');