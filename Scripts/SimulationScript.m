close all;clearvars;clc
load('Parameters.mat')%load the parameters set from the fitting
load('Data.mat');%load the real data
for i=1:length(Parameters)
par=Parameters(i,:)
alpha1  = par(1); 
beta1   = par(2);
gamma1  = par(3);
delta1  = par(4);
Lambda1 = par(5:6); 
K1      = par(7:8);
param=[alpha1,beta1,gamma1,delta1,Lambda1,K1];
dt = 0.1; % time step
% time1 = datetime(time(1)):dt:time(length(time));
time1 = datetime(time(init)):dt:datetime(2020,12,31,0,0,0);
LDdate=datefind(datetime('DD-MM-YYYY'), time1);%lock down date
N = numel(time1);
t = [0:N-1].*dt;
Init=30;Pass=15;End=90;
for conf=Init:Pass:End
% conf=1;
rel=LDdate+conf/dt;
[Y1] = SEIQRDC(param,X0,t(1:rel))
[S1;E1;I1;Q1;R1;D1;C1] =Y1; 
param1=[alpha1,beta1,gamma1,delta1,Lambda1,K1,1];
X0=Y1(end,:)
[Y2] = SEIQRDC(param1,X0,t(rel+1:end))
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
% save('Simulations.mat','S','E','I','Q','R','D','C');%for save the simulations
