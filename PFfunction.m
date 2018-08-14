% %function [ PF ] = PFfunction( t )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fc=4000 Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = 40*pi;
% P0 = 101325;  % 1 bar
% Ps = 0.5*P0;
% fc = 4000;  %Hz
% omega_c = 2*pi*fc;
% twave=2*beta/omega_c;   %there might be a *2 here
% 
% t = linspace(0,twave,10000); 
% PF= zeros(10000,1); 
% 
% 
% for i=1:10000
%    % PF(i)= 0.5*101325*exp(-4*200*200*log(10)*(t(i)-0.005)*(t(i)-0.005))*sin(8000*3.14159*(t(i)-0.005));
% 
%     PF(i)= Ps * exp(-4*(omega_c/beta)^2*log(10)*(t(i)-beta/omega_c)^2) * sin(omega_c*((t(i)-beta/omega_c)));
% end

% plot(t, PF);
%xrange(0,1.6e5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fc=1000 Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = 8*pi;
% P0 = 101325;  % 1 bar
% Ps = 0.1*P0;
% fc = 1000;  %Hz
% omega_c = 2*pi*fc;
% twave=2*beta/omega_c;   %there might be a *2 here
% 
% t = linspace(0,twave,10000); 
% PF= zeros(10000,1); 
% 
% 
% for i=1:10000
%    %PF(i)= 0.5*101325*exp(-4*250*250*log(10)*(t(i)-0.004)*(t(i)-0.004))*sin(2000*3.14159*(t(i)-0.004));
%     PF(i)= Ps * exp(-4*(omega_c/beta)^2*log(10)*(t(i)-beta/omega_c)^2) * sin(omega_c*((t(i)-beta/omega_c)));
% end
% 
% plot(t, PF/101325);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fc=3050 Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 40*pi;
P0 = 101325;  % 1 bar
Ps = 0.1*P0;
fc = 4000;  %Hz
omega_c = 2*pi*fc;
twave=2*beta/omega_c;   %there might be a *2 here
tmax=10*twave;




t = linspace(0,tmax,10000); 
PF= zeros(10000,1); 


for i=1:10000
   %PF(i)= 0.1*101325*exp(-4*200*200*log(10)*(t(i)-0.005)*(t(i)-0.005))*sin(6100*3.14159*(t(i)-0.005));
   PF(i)= Ps * exp(-4*(omega_c/beta)^2*log(10)*(t(i)-beta/omega_c)^2) * sin(omega_c*((t(i)-beta/omega_c)));
  % PF(i)=0.1*101325.0*exp(-4.0*250*250*log(10.0)*(t(i)-0.004)*(t(i)-0.004))*sin(5000*2*3.14159*(t(i)-0.004));  %5000Hz
  % PF(i)=0.1*101325.0*exp(-4.0*50*50*log(10.0)*(t(i)-0.02)*(t(i)-0.02))*sin(1000*2*3.14159*(t(i)-0.02));  %1000Hz
  % PF(i)=0.1*101325.0*exp(-4.0*5*5*log(10.0)*(t(i)-0.2)*(t(i)-0.2))*sin(100*2*3.14159*(t(i)-0.2));  %100Hz
  % PF(i)=0.1*101325.0*exp(-4.0*2.5*2.5*log(10.0)*(t(i)-0.4)*(t(i)-0.4))*sin(50*2*3.14159*(t(i)-0.4));  %50Hz
   % 
   % 
   %PF(i)=sin(6400*pi*t(i));
end
% 
figure(1)
plot(t, PF/101325);
xlim([0 0.02])
xlabel('time (s)')
ylabel('pressure (bar)')

% L = length(t);      % Signal length
% n = 2^nextpow2(L);
% 
% Fs = 10000/twave;           % Sampling frequency
% Y = fft(PF,n);
% f = Fs*(0:(n/2))/n;
% P = 2*abs(Y)/n;

NF=length(t);
FS=1/t(2);
Df=FS/NF;
f=0:Df:FS-Df;
P=2*abs( fft(PF,NF) )/NF;

figure(2)
plot(f,P./max(P)) 
%title('Gaussian Pulse in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Dimensionless amplitude')
xlim([0 8000])

