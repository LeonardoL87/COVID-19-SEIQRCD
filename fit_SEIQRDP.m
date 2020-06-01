function [alpha1,beta1,gamma1,delta1,Lambda1,K1,output,jacobian,residual] = fit_SEIQRDP(Q,R,D,Npop,E0,I0,time,guess,varargin)

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  %  option for optimset
p.addOptional('tolFun',1e-4);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;


%% Options for lsqcurvfit

% options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',1000,'Display',Display);
options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',1000);
%% Fitting the data

% Write the target input into a matrix
input = [Q;R;D];

if size(time,1)>size(time,2) && size(time,2)==1,    time = time';end
if size(time,1)>1 && size(time,2)>1,  error('Time should be a vector');end


tTarget = datenum(time-time(1)); % Number of days
t = tTarget(1):0.1:tTarget(end); % oversample to ensure that the algorithm converges
dt = median(diff(t)); % get time step


modelFun1 = @SEIQRDP_for_fitting; 

% call Lsqcurvefit
% [Coeff,resnorm,residual,exitflag,output,Lambda,jacobian]= lsqcurvefit(@(para,t) modelFun1(para,t),...
%     guess,tTarget(:)',input,zeros(1,numel(guess)),[1 2 2 2 1 1 1 1],options);
[Coeff,resnorm,residual,exitflag,output,Lambda,jacobian] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess,tTarget(:)',input,zeros(1,numel(guess)),...
    [1 2 1 2 1 2 1 2],options); 
%    [0.01 1 1/5 1/10 1 1 1 1],options);%Figuras paises (Fig 3 y A3)
%     [1 2 1 2 1 2 1 2],options);
%     [.013 2 1 2 1 2 1 2],options);
%     [0.01 1 1/5 1/10 1 1 1 1],options);
%       
%     [0.003 4 2 2 1 2 1 2],options);
    
% [Coeff] = lsqcurvefit(@(para,t) modelFun1(para,t),...
%     guess,tTarget(:)',input,zeros(1,numel(guess)),[0.025 2 1 2 1 1 1 0.09],options);
%     guess,tTarget(:)',input,zeros(1,numel(guess)),[0.003 4 2 2 1 2 1 2],options);

%% Write the fitted coeff in the outputs
alpha1 = abs(Coeff(1));
beta1 = abs(Coeff(2));
gamma1 = abs(Coeff(3));
delta1 = abs(Coeff(4));
Lambda1 = abs(Coeff(5:6));
K1 = abs(Coeff(7:8));
% if alpha1>0.07
%     alpha1=0.07;
% end
%% nested functions

    function [output] = SEIQRDP_for_fitting(para,t0)

        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5:6));
        K0 = abs(para(7:8));

        
        %% Initial conditions
        N = numel(t);
        Y = zeros(7,N);
        Y(1,1) = Npop-Q(1)-R(1)-D(1)-E0-I0;
        Y(2,1) = E0;
        Y(3,1) = I0;
        Y(4,1) = Q(1);
        Y(5,1) = R(1);
        Y(6,1) = D(1);
%         if round(sum(Y(:,1))-Npop)~=0
%             error('the sum must be zero because the total population (including the deads) is assumed constant');
%         end
        %%
        modelFun = @(Y,A,F) A*Y + F;
        dt = median(diff(t));
        
         lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
         kappa = K0(1)*exp(-K0(2).*t); % I use these functions for illustrative purpose only    
        % ODE reYution
        for ii=1:N-1
            A = getA(alpha,gamma,delta,lambda(ii),kappa(ii));
            SI = Y(1,ii)*Y(3,ii);
            F = zeros(7,1);
            F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
            Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
        end
        
%         IO = Y(3,1:N);
        Q0 = Y(4,1:N);
        R0 = Y(5,1:N);
        D0 = Y(6,1:N);
        
%         IO = interp1(t,IO,t0);
        Q0 = interp1(t,Q0,t0);
        R0 = interp1(t,R0,t0);
        D0 = interp1(t,D0,t0);
        
        output = [Q0;R0;D0];
        
        
    end
    function [A] = getA(alpha,gamma,delta,lambda,kappa)
        
        A = zeros(7);
        % S
        A(1,1) = -alpha;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-delta];
        % Q
        A(4,3:4) = [delta,-kappa-lambda];
        % R
        A(5,4) = lambda;
        % D
        A(6,4) = kappa;
        % P
        A(7,1) = alpha;
        
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end

end

