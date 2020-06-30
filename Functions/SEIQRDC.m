function [Yout] = SEIQRDC(param,X0,t)%The definition of the function was simplified for simplification propose
	alpha=param(1);
	beta=param(2);
	gamma=param(3);
	delta=param(4);
	lambda0=param(5:6);
	k0=param(7:8);
	T0=param(9);
	N = numel(t);
	Y = zeros(7,N);
	Npop=X0(7);
	Y(1,1) = Npop-X0(3)-X0(1)-X0(4)-X0(5)-X0(2)-X0(6);%originally, the protected population wasn’t taken into account for the calculation of initial Susceptible population. 
	Y(2,1) = X0(1);%E0
	Y(3,1) = X0(2);%I0
	Y(4,1) = X0(3);%Q0
	Y(5,1) = X0(4);%R0
	Y(6,1) = X0(5);%D0
	Y(7,1) = X0(6);%C0
	modelFun = @(Y,A,F) A*Y + F;
	dt = median(diff(t));
	lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); 
	k = k0(1)*exp(-k0(2).*t); T=zeros(1,N);
	alpha= alpha:-alpha/N:0-alpha/N;
	T=zeros(1,N);
	%This is the core of the model, the deconfinament process, here we define the T variable (Tau) as a time function dependent of the confinement fitted value, but also we re-define the alpha function in order to stop confinement at time T0
	if length(param)>=9 %This is done for fitting propose
	 for i=1:N
	    if i>T0/dt
	%         T(i)=1;
	%         T(i)=alpha(1)+alpha(1)*.5;
		T(i)=alpha(1)+alpha(1)*.3;
	%         T(i)=alpha(1)+alpha(1)*.1; 
		alpha(i)=alpha(i)*1/10;
	    end
	 end
	end
	for ii=1:N-1
	    V=[gamma,delta,lambda(ii),k(ii),T(ii)];
	    A = getMatrix(V,Npop);
	    Y(3,ii)=Y(3,ii);
	    SI = Y(1,ii)*Y(3,ii);
	    F = zeros(7,1);
	    F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
	    Y(:,ii+1) = RKutta(modelFun,Y(:,ii),A,F,dt,0.25,4);
	end
	S = Y(1,1:N);
	E = Y(2,1:N);
	I = Y(3,1:N);
	Q = Y(4,1:N);
	R = Y(5,1:N);
	D = Y(6,1:N);
	C = Y(7,1:N);
	[Yout]=[S;E;I;Q;R;D;C]
end

