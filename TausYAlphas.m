close all
clearvars
clc


% E0 =16; % Initial number of exposed cases. Unknown but unlikely to be zero.
% I0 = 16; % Initial number of infectious cases. Unknown but unlikely to be zero.
% Q0 = 16;
% R0 = 0;
% D0 = 0;
% P0 = 0;
I0 = 100; 
E0 =I0; % Initial number of exposed cases. Unknown but unlikely to be zero.
% I0 = 60; % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = I0;
R0 = 0;
D0 = 0;
P0 = 0;
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,output] = ...
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,output,jacobian,residual]=...
%     fit_SEIQRDP(Confirmed(init:limit)-Recovered(init:limit)-Deaths(init:limit),...
%     Recovered(init:limit),Deaths(init:limit),Npop,E0,I0,time(init:limit),guess);

% %fited parameters (Spain)
% % ---------------------
% alpha1  = 0.012999999999976  ;
% beta1   = 0.999999734304382  ;
% gamma1  = 0.999999983687989  ;
% delta1  = 0.532570498284943  ;
% Lambda1 =[0.0650722743882673 0.0386318696954017];
% Kappa1  =[0.0129036440278137 1.13687095960597e-09];
alpha1  = 0.012;
beta1   = 1.1027;
gamma1  = 0.9188;
delta1  = 0.5301;
Lambda1 =[ 0.0507,0.0243];
Kappa1  =[0.0687,0.0370];
% % ----------------------

Npop=50e6;
dt = 1; % time step

time1 = datetime(datetime(2020,2,1,0,0,0):dt:datetime(2021,12,31,0,0,0));
fecha=datefind(datetime('1-Mar-2020'), time1);

N = numel(time1);
t = [0:N-1].*dt;

Init=30;Pass=15;End=60;
k=1;

ini1=10;inc1=30;fin1=200;
ini2=0.1;inc2=0.2;fin2=.9;
c=1;
% alpha1=alpha1*2;   
for i=ini2: inc2: fin2 
    ratio=i;
    fprintf('Iter: %i de %i, ratio %i\n',k,length(ini2: inc2: fin2),ratio)
    tau=alpha1/ratio;
    Soltados=0;
    for conf=Init:Pass:End  
        suelto=fecha+(conf)/dt;
        % conf=1;
        [S1,E1,I1,Q1,R1,D1,P1] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,...
                    E0,I0,Q0,R0,D0,P0,t(1:suelto));
        [S2,E2,I2,Q2,R2,D2,P2] = SEIQRDP3(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,conf,tau,...
                E1(end),I1(end),Q1(end),R1(end),D1(end),P1(end),t(suelto+1:end));     
        S=[S1,S2];E=[E1,E2];I=[I1,I2];Q=[Q1,Q2];R=[R1,R2]; D=[D1,D2];P=[P1,P2];
%         fprintf('Individos liberados %i\n',D2(end)-D1(end))
%         Soltados=Soltados+P1(end)-P2(end);
        switch c
                case 1
        %             conf=30;
                    color='k';
%                     Soltados=Soltados+P1(end)-P2(end);
                case 2
        %             Subfolder='45Days/';
        %             conf=45;
                    color='b';
%                      Soltados=Soltados+P1(end)-P2(end);
                case 3
        %             Subfolder='60Days/';
        %             conf=60;
                    color='r';
                     Soltados=Soltados+P1(end)-P2(end);
        end
        
        figure(k)
        subplot(2,1,1)
        plot(time1,Q+I,color,'LineWidth',1)
        hold on
        ylabel('Cases')
%             str = sprintf('Infected $\\frac{\\alpha}{\\tau}=%.2f$. $\\alpha= %.3f \\tau=%.3f$',alpha1/tau,alpha1,tau);
        str = sprintf('Infected $\\frac{\\alpha}{\\tau}=%.2f$',alpha1/tau);

        title(str,'Interpreter','latex')
        axis tight
        set(gca,'yscale','lin')
        str = sprintf('days = %i',conf);
        [val,pos]=max(Q);
        text(time1(pos),Q(pos),str, 'FontSize',7)

        grid minor


        Cases=Q+D+R+I;
        subplot(2,1,2)
        plot(time1,Cases,color,'LineWidth',2)
        hold on
        ylabel('Total Cases')
        xlabel('Time (days)')
        set(gcf,'color','w')
        axis tight
        set(gca,'yscale','lin')
        str = sprintf('days = %i',conf);
        text(time1(end),Cases(end),str, 'FontSize',7);
        str = sprintf('Total Infected $\\frac{\\alpha}{\\tau}=%.2f$',alpha1/tau);
        title(str,'Interpreter','latex')
        grid minor
        c=c+1;

    end
     fprintf('Individos liberados %i\n',Soltados)
        c=1;
%         str = sprintf('Ratio=%.2f.eps',alpha1/tau);
%         saveas(gcf,str,'epsc')
    k=k+1;
%         close all
end

% Tau=1./(1+exp(-t))*alpha1*1.3;
