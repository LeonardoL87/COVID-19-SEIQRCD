function [parout,output,jacobian,residual] = fit_Model(X0,X1,time,guess)
	Q=X0(1,:);
	R=X0(2,:);
	D=X0(3,:);
	Npop=X1(1);
	E0=X1(2);
	I0=X1(3);
	p = inputParser();
	p.CaseSensitive = false;
	p.addOptional('tolX',1e-12);  %  less tolerance for 
	p.addOptional('tolFun',1e-12);  %  option for optimset
	tolX = p.Results.tolX ;
	tolFun = p.Results.tolFun ;
	% Options for the optimizer
	options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',10000);%more evaluations in order to not end the fitting before reach the minimum
	% Write the target input into a matrix
	IN = [Q;R;D];
	tTarget = datenum(time-time(1)); 
	t = tTarget(1):0.1:tTarget(end);
	dt = median(diff(t)); 
	up = [1 1 1 1 1 1 1 1];
	modelFun1 = @SEIQRDCFitting; 
	[Coeff,resnorm,residual,exitflag,output,Lambda,jacobian] =lsqcurvefit(@(para,t) modelFun1(para,t),guess,tTarget(:)',IN,zeros(1,numel(guess)),ub,options);
	% Write the fitted parameters
	parout(1) = abs(Coeff(1));
	parout(2) = abs(Coeff(2));
	parout(3) = abs(Coeff(3));
	parout(4) = abs(Coeff(4));
	parout(5:6) = abs(Coeff(5:6));
	parout(7:8) = abs(Coeff(7:8));

	function [output] = SEIQRDCFitting(para,t0)%for fitting propose this functions calls the main function and return the populations of interest for the Error calculation
		X=zeros(1,7);
		X(1)=E0;
		X(2)=I0;
		X(3)=Q(end);
		X(4)=R(end);
		X(5)=D(end);
		X(6)=0;
		X(7)=Npop;
		[Yout] = SEIQRDC(para,X,t0)
		Q0 = interp1(t,Yout(4,1:N),t0);
		R0 = interp1(t,Yout(5,1:N),t0);
		D0 = interp1(t,Yout(6,1:N),t0);
		output = [Q0;R0;D0];
	end
end

