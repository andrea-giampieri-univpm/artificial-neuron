clear all
close all
clc

n = 1000;
T = 0.01;
h=0.01;
tempo = zeros(1,n);
Vc = zeros(1,n);
Ana = zeros(1,n);
Ak = zeros(1,n);
Ip = zeros(1,n);
Ikp = zeros(1,n);
Vnap = zeros(1,n); 
Inap = zeros(1,n);
Ia = zeros(1,n);
Is = zeros(1,n);
Ic = zeros(1,n);
Ika = zeros(1,n);
Inaa = zeros(1,n);
Iext = zeros(1,n);

%Dati sul circuito: 
Gna = 0.17;
Dna = -0.06;
i1 = 0.1;
i2 = 0.3;
Ena = 0.6;
Gk = 1;
Dk = -1.25;
v1 = 0.5;
v2 = 2;
Ek = -0.7;
C = 0.01;
l = 0.05;
y = 0.1; Ga=1/y;
e = 0.01;

%parametri soglia
Imthr=(Gna+Gk+Ga)*(Ek+v1)-(Gna*Ena+Gk*Ek);
Inathr=(i1*(1+y*Gna)/(y*Gna))+(Ena/y);
Ikthr=(v1*y/(1+Gk))+(Ena/y);

im=(i1+i2)/2;%corrente media canale sodio
p=sqrt(abs(Dna)/(Gna+Dna))*(i2-i1)/2; %parametro rho canale sodio
vm=(v1+v2)/2; %tensione media canale potassio
u=sqrt(Gk/abs(Gk+Dk))*(v2-v1)/2; %parametro mu canale potassio
 
syms Vna(I) Ik(V)
 
%tensione canale sodio S-shape
%Vna(I) = I/Gna + p/Dna * (((atan((I-im)/p))) + ((atan(im/p))));
Vna(I) = piecewise((I<=i1),I/Gna,(i1<I) & (I<i2),I/Gna +(I-i1)/Dna,(I>=i2),I/Gna + (i2-i1)/Dna);

%corrente canale potassio N-shape
%Ik(V) = Gk*V+Dk*u*((atan((V-vm)/u))) + (Dk*u*((atan(vm/u))));
Ik(V) = piecewise((V<=v1),Gk*V,(v1<V) & (V<v2),Gk*V +(V-v1)*Dk,V>=v2,Gk*V +(v2-v1)*Dk );
%corrente simulazione
Iext([1:800])=0;
Iext([801:1200])=7;
Iext([1201:2000])=0;

%valori iniziali
% Vc(1)=Ena*(Gna/(Gna+Ga))+Iext(1)/(Gna+Ga);
%  Ana(1)=Ga*Vc(1);
% Ekv=Gk*Ek/(Gk+Ga)+Iext(1)/(Gk+Ga);
% Ak(1)=-Ga*Ekv;
% Inap(1)=-Ana(1)-Ak(1)+Iext(1);
% Ikp(1)=-Ga*Ek+Iext(1);
% Ikp(1)=Ik(Vc(1)-Ek);

Ak(1)=3; %la cond iniziale aumenta la frequenza delle spike
Ana(1)=0.1; 
Inaa(1)=Inap(1)+Ana(1); %somma dei canali sodio
Ika(1)=Ikp(1)+Ak(1); %somma dei canali potassio
Ia(1)=Ana(1)-Ak(1); %corrente attiva pompa ionica
Is(1)=Ana(1)+Ak(1); %somma correnti pompa ionica

for i=1:n-1
    100*i/n %percentuale
    tempo(i+1)=tempo(i)+T;

    %sistema
    vF1=(-Inap(i)-Ikp(i)-Ana(i)+Ak(i)+Iext(i))/C;
    vF2=(-(Inap(i)+h*vF1)-(Ikp(i)+h*vF1)-(Ana(i)+h*vF1)+(Ak(i)+h*vF1)+(Iext(i)+h*vF1))/C;
    Vc(i+1) = Vc(i)+(vF1)*h/2+(vF2)*h/2;
    
    naF1=l*Ana(i)*(Vc(i)-y*(Ana(i)-Ak(i)));
    naF2=l*(Ana(i)+h*naF1)*((Vc(i)+h*naF1)-y*((Ana(i)+h*naF1)-(Ak(i)+h*naF1)));
    Ana(i+1)=Ana(i)+(naF1)*h/2+(naF2)*h/2;
    
    kF1=l*Ak(i)*(-Vc(i)+y*(Ana(i)-Ak(i)));
    kF2=l*(Ak(i)+h*kF1)*(-(Vc(i)+h*kF1)+y*((Ana(i)+h*kF1)-(Ak(i)+h*kF1)));
    Ak(i+1)=Ak(i)+(kF1)*h/2+(kF2)*h/2;

    Ikp(i+1)=Ik(Vc(i+1)-Ek);
    Inap(i+1)=Gna*(Vc(i+1)-Ena);
    Vnap(i+1)=Vna(Inap(i+1));

    %valori plot
    Inaa(i+1)=Inap(1+i)+Ana(i+1); %somma dei canali sodio
    Ika(i+1)=Ikp(1+i)+Ak(i+1); %somma dei canali potassio
    Ia(i+1)=Ana(i+1)-Ak(i+1); %corrente attiva pompa ionica
    Is(i+1)=Ana(i+1)+Ak(i+1); %somma correnti pompa ionica
end

%correnti
Inap(n)=Gna*(Vc(n)-Ena);
Vnap(n)=Vna(Inap(n));    
Ikp(n)=Ik(Vc(n)-Ek);

subplot(3,4,1);
fplot (Ik); hold on;
title({'','Ik / Vk',''});
subplot(3,4,2);
fplot (Vna)
title({'','Vna / Ina',''});
subplot(3,4,3);
plot (tempo,Ikp)
title({'','Ikp / tempo',''});
subplot(3,4,4);
plot (tempo,Inap)
title({'','Inap / tempo',''});
subplot(3,4,5);
plot (tempo,Vnap)
title({'','Vnap / tempo',''});
subplot(3,4,6);
plot (tempo,Ana)
title({'','Ana / tempo',''});
subplot(3,4,7);
plot (tempo,Ak)
title({'','Ak / tempo',''});
subplot(3,4,8);
plot (tempo,Vc)
title({'','Vc / tempo',''});
subplot(3,4,9);
plot (tempo,Ia)
title({'','Ia / tempo',''});
subplot(3,4,10);
plot (tempo,Is)
title({'','Is / tempo',''});
subplot(3,4,11);
plot (tempo,Ika)
title({'','Ik all / tempo',''});
subplot(3,4,12);
plot (tempo,Inaa)
title({'','INa all / tempo',''});