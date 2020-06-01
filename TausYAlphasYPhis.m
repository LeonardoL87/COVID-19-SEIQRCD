close all
clearvars
clc


I0 = 1e4; % Initial number of infectious cases. Unknown but unlikely to be zero.
E0 =I0; % Initial number of exposed cases. Unknown but unlikely to be zero.
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
alpha1  = 0.015;
beta1   = 1.2027;
gamma1  = 0.9188;
delta1  = 0.5301;
Lambda1 =[ 0.0507,0.0243];
Kappa1  =[0.0687,0.0370];
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,output] = ...
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,output,jacobian,residual]=...
%     fit_SEIQRDP(Confirmed(init:limit)-Recovered(init:limit)-Deaths(init:limit),...
%     Recovered(init:limit),Deaths(init:limit),Npop,E0,I0,time(init:limit),guess);

% %fited parameters (Spain)
% % ---------------------
% alpha1  = 0.01;
% beta1   = 1.9994;
% gamma1  = 0.9989;
% delta1  = 1.1088;
% Lambda1 =[ 0.9689,0.0016];
% Kappa1  =[0.0162,0.0001];
% % ----------------------

Npop=50e6;
dt = 1; % time step

time1 = datetime(datetime(2020,2,1,0,0,0):dt:datetime(2021,12,31,0,0,0));
fecha=datefind(datetime('1-Mar-2020'), time1);


N = numel(time1);
t = [0:N-1].*dt;

Init=30;Pass=15;End=60;
k=1;

% ini1=30;inc1=30;fin1=120;
ini2=120;inc2=60;fin2=365;
% for j=ini1: inc1: fin1
%     alpha1=j;
    
c=1;      
for i=ini2: inc2: fin2
        phi=1/i;
        fprintf('Iter: %i de %i\n',k,length(ini2: inc2: fin2))
        ratio=.5;
        tau=alpha1/ratio;
%         tau=alpha1+alpha1*.1;
        for conf=Init:Pass:End
            suelto=fecha+(conf)/dt;
            [S1,E1,I1,Q1,R1,D1,P1] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,...
                    E0,I0,Q0,R0,D0,P0,t(1:suelto));
                
            [S2,E2,I2,Q2,R2,D2,P2] = SEIQRDP4(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,conf,tau,phi,...
                    E1(end),I1(end),Q1(end),R1(end),D1(end),P1(end),t(suelto+1:end)); 
    %         [S,E,I,Q,R,D,P] = SEIQRDP5(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,conf,tau,i,...
    %                 E0(end),I0(end),Q0(end),R0(end),D0(end),P0(end),t);

            S=[S1,S2];E=[E1,E2];I=[I1,I2];Q=[Q1,Q2];R=[R1,R2]; D=[D1,D2];P=[P1,P2];
            % conf=1;
%             [S,E,I,Q,R,D,P] = SEIQRDP3(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,conf,tau,E0,I0,Q0,R0,D0,P0,t);
%             [S,E,I,Q,R,D,P] = SEIQRDP4(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,conf,tau,phi,E0,I0,Q0,R0,D0,P0,t);
            
            switch c
                case 1
        %             conf=30;
                    color='k';
                case 2
        %             Subfolder='45Days/';
        %             conf=45;
                    color='b';
                case 3
        %             Subfolder='60Days/';
        %             conf=60;
                    color='r';
            end

            figure(k)
            subplot(2,1,1)
            plot(time1,Q+I,color,'LineWidth',2)
            hold on
            ylabel('Cases')
%             str = sprintf('Infected $\\frac{\\alpha}{\\tau}=%.2f$. $\\alpha= %.3f \\tau=%.3f$',alpha1/tau,alpha1,tau);
            str = sprintf('Infected $\\frac{\\alpha}{\\tau}=%.2f$. $\\phi=%i$',ratio,i);
            
            title(str,'Interpreter','latex')
            axis tight
            set(gca,'yscale','lin')
            str = sprintf('days = %i',conf);
            [val,pos]=max(Q);
            text(time1(pos),Q(pos),str, 'FontSize',7)

            grid minor
            
           for l=2:length(R)
                if R(l)<R(l-1)
                    R(l)=R(l-1);
                end
            end


            Cases=Q+D+R;
%             Cases=Acum(Q)*1e-2;%-Acum(R)-Acum(D);
%             Cases=Acum(Q)-Acum(R)-Acum(D);
            for l=2:length(Cases)
                if Cases(l)<Cases(l-1)
                    Cases(l)=Cases(l-1);
                end
            end
            [val,loc]=max(Cases);
            Cases(loc:end)=val;
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
            str = sprintf('Total Infected $\\frac{\\alpha}{\\tau}=%.2f$. $\\phi=%i$',ratio,i);
            title(str,'Interpreter','latex')
            grid minor
            
            c=c+1;

        end
        c=1;
        str = sprintf('Ratio=%.2finmunidad=%idias.eps',ratio,i);
        saveas(gcf,str,'epsc')
        k=k+1;
%         close all
end
% end

% Tau=1./(1+exp(-t))*alpha1*1.3;
