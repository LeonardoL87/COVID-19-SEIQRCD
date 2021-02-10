clearvars;close all;clc;
addpath('The Path')
load('Data.mat');%Data it's the real data fro the location you want to fit
init=1;
%In case you want to select an specific date rage
% limit=datefind(datetime('13-Mar-2020'), time);
% init=datefind(datetime('13-Mar-2020'), time);
limit=length(time);
iter=1000;%number of iterations
Paramters=zeros(iter,12);
for i=1:iter
% Fitting of the generalized SEIR model to the real data
% Definition of the first estimates for the parameters
alpha0 =(.09-.01).*rand(1,2) + .01;% random number between 0.01 and 0.09; %protection rate
beta0= rand(1,2);%1.0; % Infection rate
LT0 = (10-1).*rand() + 1;% random number between 1 and 10; % latent time in days
QT0 = (30-14).*rand() + 14;% random number between 14 and 30; % quarantine time in days
lambda0 = rand(1,2);%[0.1,0.05]; % recovery rate
K0 = rand(1,2);%[0.1,0.05]; % death rate
guess = [alpha0,beta0,1/LT0,1/QT0,lambda0,K0];
% Set Initial conditions: All the vectors are contained in Data.mat
E0 = Confirmed(init); % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 =  Confirmed(init); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = Confirmed(init);
R0 = Recovered(init);
D0 = Deaths(init);
P0 = 0;
QD=Confirmed(init:limit)-Recovered(init:limit)-Deaths(init:limit);
RD=Recovered(init:limit);
DD=Deaths(init:limit);
X0=[QD;RD;DD];
X1=[Npop;E0;I0]
fprintf('Optimzation N°: %i of %i\n',i,iter)
[parout,output,jacobian,residual] = fit_Model(X0,X1,time(init:limit),guess)
%Show the output of the process    
output
Paramters(i,:)=[parout(1),parout(2),parout(3),parout(4),parout(5:6),parout(7:8)];    
end
% save('Paramters_optim.mat','Paramters');%for save the parameters
