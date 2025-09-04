% Calculate responses of the benchamek building with rayleigh damping 
clear all;clc

dt=input('Time step  ');

% Input forces
%S=input('Force intensity of white noise   ');
S=20.0;

%Duration=input('Time duration  ');
Duration=5;
SeedNum=input('Input the seed number  ');
%Filter index
Findx=1;
%Filter index =1 -  use a filter to handle the direct pass-through problem. 

%Noise level
nl=input('noise level  ');

% ***** Generate white noise
randn('state',SeedNum);
t=[dt:dt:Duration]';
Nt=length(t);

if Findx==1
   filter_order = 6;
   filter_cutoff = 100; %Hz
   [filt_num,filt_den] = butter(filter_order,filter_cutoff*2*dt);
   Nt2=Nt+2*filter_order;
else
   Nt2=Nt;
end;

ff=S*randn(Nt2,1)./dt^0.5;

if Findx==1
   ff= filter(filt_num,filt_den,ff);
   ff= ff(Nt2-Nt+1:end,:);
end;

force(1,:)=0;
force(2:Nt+1,1)=t;
force(2:Nt+1,2)=ff;

save a4.txt force -ascii