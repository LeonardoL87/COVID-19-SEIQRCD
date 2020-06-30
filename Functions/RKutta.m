function [Y] = RKutta(Fun,Y,A,F,dt,h,O) %The integration function was improved in order to incorporate a variable size-step h, which is important for reproduce some dynamics 
      if O==4
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*h*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*h*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*h*dt,A,F);
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
      	else % Runge-Kutta of order 2
	if O==2
	   k_1 = Fun(Y,A,F);
   	   k_2 = Fun(Y+0.5*h*dt*k_1,A,F);
   	   Y = Y + (1/6)*(k_1+2*k_2+)*dt;
	end
      end
end

