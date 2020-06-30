close all;clearvars;clc
load('Parameters.mat')
load('Data.mat')
%Set initial conditions
I0 = 100; 
E0 =I0; 
Q0 = I0;
R0 = 0;
D0 = 0;
P0 = 0;
par=mean(Parameters);
alpha1  = par(1);
beta1   = par(2);
gamma1  = par(3);
delta1  = par(4);
Lambda1 =par(5:6);
Kappa1  = par(7:8);
dt =0.1; % time step
time1 = datetime(datetime(2020,2,1,0,0,0):dt:datetime(2021,12,31,0,0,0));
LDdate=datefind(datetime('DD-MM-YYYY'), time1);
N = numel(time1);
t = [0:N-1].*dt;
Init=30;Pass=15;End=60;
k=1;
ini2=120;inc2=60;fin2=356;
c=1;
for i=ini2: inc2: fin2 
ratio=0.5;
fprintf('Iter: %i de %i\n',k,length(ini2: inc2: fin2))
tau=alpha1/ratio;
for conf=Init:Pass:End
rel=LDdate+(conf)/dt;
param=[alpha1,beta1,gamma1,delta1,Lambda1,K1];
[Y1] = SEIQRDC(param,X0,t(1:rel));
[S1;E1;I1;Q1;R1;D1;C1] =Y1; 
X0=Y1(end,:)
param1=[alpha1,beta1,gamma1,delta1,Lambda1,K1,1,tau,i/dt];
[Y2] = SEIQRDCPhi(param1,X0,t(rel+1:end));%this is a version of SEIQRDC in which we can set the value of tau according to the ratio. The only difference with the reported in this document it's T(i)= tau and the parameter phi is introduced in the system (i). In this way, the GetMatrix function incorporates  the entries that model phi=i;
[S2;E2;I2;Q2;R2;D2;C2] =Y2; 
S=[S1,S2];
E=[E1,E2];
I=[I1,I2];
Q=[Q1,Q2];
R=[R1,R2];
D=[D1,D2];
C=[C1,C2];
end
end
% save('SimulationsImmunityLoss.mat','S','E','I','Q','R','D','C');%To save the simulations  


