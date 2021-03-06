function [M] = getMatrix(Coeff,Npop)
        alpha=Coeff(1);
        gamma=Coeff(2);
        delta=Coeff(3);
        lambda=Coeff(4);
        k=Coeff(5); 
        T=Coeff(6);
        phi=0;%if no loss of immunuty is simulated then re-susceptibility term becomes 0
        mu=1/(80*365);
% For temporalimmunity
  if length(Coeff>6)
      phi=Coeff(7);
  end 
syms S E I Q R D C 
vars=[S E I Q R D C];
eqns = [mu*Npop + tau*C - alpha * S-mu*S + phi*R;
        gamma*E-mu*E;
        gamma*E- delta*I-mu*I;
        delta*I-lambda*Q-k*Q-mu*Q;
        lambda*Q-mu*R- phi*R;
        k*Q;
        alpha*S-mu*C-tau*C];
[M,b] = equationsToMatrix(eqns,vars);
end

