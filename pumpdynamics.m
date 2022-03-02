

close all
clear all
clc
% Array per la simulazione e campionamento

n = 1000; % Numero di campioni temporali
T = 0.01; % Intervallo di campionamento
h= 0.01;
tempo=zeros(1,n);
Ic=zeros(1,n); %corrente condensatore
Ip=zeros(1,n); %corrente pompa
Ana=zeros(1,n); %corrente pompa sodio
AnaH=zeros(1,n); %corrente pompa sodio

Ak=zeros(1,n);
AkH=zeros(1,n);

IcH=zeros(1,n); %corrente condensatore


l= 0.1; %parametro lambda
y=0.1; %parametro gamma
C=0.01;%condensatore
Vc=[-5:0.01:5]; %vettore tensioni
%Vc=zeros(1,n);

Ana(1)=0.5;
AnaH(1)=0.5;
Ak(1)=0.1;
AkH(1)=0.1;
 for i=1:n
    tempo(i+1)=tempo(i)+T;
    Ic(i)=Ak(i)-Ana(i);
    IcH(i)=AkH(i)-AnaH(i);
    
    Ana(i+1)=Ana(i)*(T*Vc(i)*l - (Ana(i)-Ak(i))*T*l*y + 1);
    Ak(i+1)=-Ak(i)*(T*Vc(i)*l + (Ana(i)-Ak(i))*T*l*y - 1);
    
    naF1=l*AnaH(i)*(Vc(i)-y*(AnaH(i)-AkH(i)));
    naF2=l*(AnaH(i)+h*naF1)*((Vc(i)+h*naF1)-y*((AnaH(i)+h*naF1)-(AkH(i)+h*naF1)));
    AnaH(i+1)=AnaH(i)+(naF1)*h/2+(naF2)*h/2;
    
    kF1=l*AkH(i)*(-Vc(i)+y*(AnaH(i)-AkH(i)));
    kF2=l*(AkH(i)+h*kF1)*(-(Vc(i)+h*kF1)+y*((AnaH(i)+h*kF1)-(AkH(i)+h*kF1)));
    AkH(i+1)=AkH(i)+(kF1)*h/2+(kF2)*h/2;
 end
Ic(n+1)=Ak(n+1)-Ana(n+1);
IcH(n+1)=AkH(n+1)-AnaH(n+1);
subplot(1,3,1);
plot (Vc,Ana,'red'); hold  on;
plot (Vc,AnaH,'blue');
title({'','Ana',''});

subplot(1,3,2);
plot (Vc,Ak,'red');hold on;
plot (Vc,AkH,'blue');
title({'','Ak',''});

subplot(1,3,3);
plot (Vc,Ic,'red');hold on;
plot (Vc,IcH,'blue');
title({'','Ic',''});