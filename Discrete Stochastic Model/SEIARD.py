#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:27:20 2020

@author: leonardo
"""
"""
Stochastic model COVID-19
stochastic version of the NHB and RINP compartmental model
Odds are deducted from rates
Populations:
     - S: Susceptible
     - E: Exposed
     - I: Infected
     - R: Recovered
     - D: Dead
     - Y: Total cases
Parameters:
      * alpha: confination rate
      * beta: infection rate
      * gamma: incubation rate
      * delta: detected infected (Not included here)
      * sigma: sympthomatic rate
      * Lambda: recovery rate
      * kappa: death rate
      * tau: deconfination rate
"""
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit, Model
import math as m
from scipy import interpolate
from sklearn.metrics import mean_squared_error 


np.random.seed(123)
#

#------------------------------------------------------------------------------
#                       Diference equation system
#               This can be formulated from the ODE system
#               Return a column vector with the system at time t
#------------------------------------------------------------------------------
def SEAIRDP(t, y, ps):   
#    alpha  = ps['alpha0'].value # for constant alpha
    alpha  = ps['alpha0'].value*m.exp(-ps['alpha1'].value*t)
    bet    = ps['beta'].value
    gamma  = ps['gamma'].value
    Lambda = ps['lambda0'].value
#    Lambda = ps['lambda0'].value*(1.-m.exp(-ps['lambda0'].value*t))
    kappa  = ps['kappa0'].value
#    kappa = ps['kappa0'].value*m.exp(-ps['kappa0'].value*t)    
#    tau    = ps['tau0'].value # for constant tau
    tau    = ps['tau0'].value*(1.-m.exp(-ps['tau0'].value*t))
#        delta = param[8]
    sigma  = ps['sigma'].value
    S,E,A,I,R,D,P,Y=y
    beta=(bet*(I + 0.5*A))/N
    
    efrac =   (1.0 - m.exp(-beta*dt))        #|exposedexposed prob
    ifrac =   (1.0 - m.exp(-gamma*dt))       #|infection prob
    rfrac =   (1.0 - m.exp(-Lambda*dt))      #|recov prob
    pfrac =   (1.0 - m.exp(-alpha*dt))       #|protec prob
    dfrac =   (1.0 - m.exp(-kappa*dt))       #|death prob
    relfrac = (1.0 - m.exp(-tau*dt))         #|release prob

    
#    --------------------------------------------

    Su_E_P =  np.random.multinomial(S,[1-pfrac,pfrac])
    exposed =  np.random.binomial(Su_E_P[0],efrac)
    protected = Su_E_P[1]
    
    Ex_I_A = np.random.multinomial(E,[1-sigma,sigma])
    infection = np.random.binomial(Ex_I_A[0],ifrac)
    asympthomatic = np.random.binomial(Ex_I_A[1],ifrac)
    
    I_R_D=np.random.multinomial(I,[1-dfrac,dfrac])

    recovery = np.random.binomial(I_R_D[0],rfrac)
    deaths = I_R_D[1]
    recoveryA = np.random.binomial(A,rfrac)
    released = np.random.binomial(P,relfrac)
#    ------------------------------------------
    S = S - exposed - protected + released                  #| Susceptible
    E = E + exposed - infection- asympthomatic              #| Exposed
    A = A + asympthomatic - recoveryA                       #| Asympthomatic
    I = I + infection - recovery - deaths                   #| Infected 
    R = R + recovery + recoveryA                            #| Recovered
    D = D + deaths                                          #| Deaths
    P = P + protected - released                            #| Protected
    Y = Y + infection 			                      #| Total Cases
#    ------------------------------------------
    return [S,E,A,I,R,D,P,Y]
                               
#------------------------------------------------------------------------------
#                       Simulation function
#               Return a matrix with the dynamic of each population
#------------------------------------------------------------------------------
def simulate(t,u,ps): 
    S = np.zeros(len(t))
    E = np.zeros(len(t))
    A = np.zeros(len(t))
    I = np.zeros(len(t))
    R = np.zeros(len(t))
    D = np.zeros(len(t))
    P = np.zeros(len(t))
    Y = np.zeros(len(t))
    for j in range(len(t)):
        u = SEAIRDP(t[j],u,ps)
        S[j],E[j],A[j],I[j],R[j],D[j],P[j],Y[j] = u
    return  {'t':t,'S':S,'E':E,'A':A,'I':I,'R':R,'D':D,'P':P,'Y':Y}
#------------------------------------------------------------------------------
#                       Simulation function N times
#               Return a matrix with the mean dynamic of each population 
#------------------------------------------------------------------------------
def  simulate_N_times(t, u, ps):  
    times=100
    S=np.zeros([times,len(t)])
    E=np.zeros([times,len(t)])
    A=np.zeros([times,len(t)])
    I=np.zeros([times,len(t)])
    R=np.zeros([times,len(t)])
    D=np.zeros([times,len(t)])
    P=np.zeros([times,len(t)])
    Z=np.zeros([times,len(t)])
    vec2=np.zeros(len(t))
#    Yi={'t':t,'S':vec,'E':vec,'A':vec,'I':vec,'Q':vec,'R':vec,'D':vec,'P':vec,'Y':vec}
    Y={'t':t,'S':vec2,'E':vec2,'A':vec2,'I':vec2,'R':vec2,'D':vec2,'P':vec2,'Y':vec2}
    
    for i in np.arange(times):
        y = simulate(t, u, ps)
        S[i,:]=y['S']
        E[i,:]=y['E']
        A[i,:]=y['A']
        I[i,:]=y['I']
#        Q[i,:]=y['Q']
        R[i,:]=y['R']
        D[i,:]=y['D']
        P[i,:]=y['P']
        Z[i,:]=y['Y']

    Y['S']=S.mean(0)
    Y['E']=E.mean(0)
    Y['A']=A.mean(0)
    Y['I']=I.mean(0)
#    Y['Q']=Q.mean(0)*times
    Y['R']=R.mean(0)
    Y['D']=D.mean(0)
    Y['P']=P.mean(0)
    Y['Y']=Z.mean(0)
#    
    return Y


#------------------------------------------------------------------------------
#                       Interpolation function
#              Single interpolation of 1 vector
#               
#------------------------------------------------------------------------------

def interpolation(y,t,ti):
    f= interpolate.interp1d(t,y, kind='nearest')
    f2 =f(ti)
    return f2
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#                       Interpolation function
#              Single interpolation of 1 vector
#               
#------------------------------------------------------------------------------

def Sys_interpolation(Y,t,ti):
    col=Y.columns
    datcol=col[1:len(col)]
    Yinterp={}
    Yinterp['t']=ti
    for i in datcol:
        yi=Y[str(i)].to_numpy()
        f2=interpolation(yi,t,ti)
        Yinterp[str(i)]=f2
    return Yinterp
#------------------------------------------------------------------------------
    

#------------------------------------------------------------------------------
#                       Residuals function
#              Get the resilduals calculating the function output for a
#               set of parameters given by the optimization method
#------------------------------------------------------------------------------

def residual(ps,ts,data):
    model0 = pd.DataFrame(simulate_N_times(ts, y0, ps))
    ti= np.linspace(ts[0],ts[len(ts)-1],int(ts[len(ts)-1]-ts[0]+1))
    model = pd.DataFrame(Sys_interpolation(model0,ts,ti))
    Q=model.I
    R=model.R
    D=model.D
    Inf=((Q - data[0,:])/Npop).ravel()
    Rec=((R - data[1,:])/Npop).ravel()
    Dea=((D - data[2,:])/Npop).ravel()
    resid=np.array([Inf,Rec,Dea])
    return resid

#===============================LOAD DATA FROM FILES===========================
#Loading Data
nombredearchivo_csv='Your Path/your file.csv'
Npop=66.99e6
Data = pd.read_csv(open(nombredearchivo_csv), sep=",")


time=pd.to_datetime(Data.Time)
ini=0
fin=150
#ini=150
#fin=len(Data)
#Active = Data.Active[ini0:fin]
#pos=np.where(Active>= 1)[0]#40
#ini=pos[0]
time = pd.to_datetime(Data.Time[ini:fin])
Recovered = Data.Recovered[ini:fin]
Deaths = Data.Deaths[ini:fin]
Confirmed = Data.Confirmed[ini:fin]
Active = Data.Active[ini:fin]

DataM=np.array([Active,Recovered,Deaths])


I0 = Active[ini]+200
E0 = I0
Q0 = I0
A0 = I0
R0 = Recovered[ini]
D0 = Deaths[ini]
P0 = 0
S0 = Npop-E0-I0-Q0-R0-D0-P0-A0
Y0 = I0+A0

#===============SETING THE ORIGIAL SET OF PARAMETERS===========================
lb=0
ub=1
params = Parameters()
params.add('alpha0', value=np.random.uniform(1e-3,1), min=lb, max=1)
params.add('alpha1', value=np.random.uniform(1e-1,1), min=lb, max=1)
params.add('beta', value= np.random.uniform(.8,2), min=lb, max=3)
params.add('gamma', value= np.random.uniform(1e-1,1), min=lb, max=1)
params.add('lambda0', value= np.random.uniform(0,0.03), min=0, max=1)
#params.add('lambda1', value= np.random.uniform(0,0.1), min=0, max=1)
params.add('kappa0', value=np.random.uniform(0,0.03), min=0, max=1)
#params.add('kappa1', value= np.random.uniform(0,0.1), min=0, max=1)
params.add('tau0', value= np.random.uniform(1e-3,1), min=lb, max=1)
params.add('tau1', value= np.random.uniform(1e-1,1), min=lb, max=ub)
params.add('sigma', value= np.random.uniform(0.6,0.7), min=lb, max=1)


N=Npop
dt=1/24

tf = len(time)
tl = int(tf/dt)
t = np.linspace(0,tf-1,tl)

y0 = [S0,E0,A0,I0,R0,D0,P0,Y0]

#===============GETING THE SOLUTION OF THE MODEL===============================

sol = minimize(residual, params, args=(t,DataM),method='least_squares',max_nfev=1000,
                   ftol=1e-12,gtol=1e-12,xtol=1e-12,loss='soft_l1',diff_step=1e-3,verbose=2,tr_solver='lsmr')

paropt=sol.params
print(report_fit(sol))
#sir_out = pd.DataFrame(simulate(paropt,t))

sir_out = pd.DataFrame(simulate_N_times(t, y0, paropt))


#sir_out = pd.DataFrame(simulate(params,t))
ti= np.linspace(t[0],t[len(t)-1],int(t[len(t)-1]-t[0])+1)
sir_out = pd.DataFrame(Sys_interpolation(sir_out,t,ti))

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(time,"I", 'r',data=sir_out, alpha=0.5, lw=2, label='Active')
ax.plot(time,"R", 'g',data=sir_out, alpha=0.5, lw=2, label='Recovered')
ax.plot(time,"D", 'c',data=sir_out, alpha=0.5, lw=2, label='Deaths')
ax.plot(time,Active, 'or', alpha=0.5, lw=2, label='Active')
ax.plot(time,Recovered, 'og', alpha=0.5, lw=2, label='Recovered')
ax.plot(time,Deaths, 'oc', alpha=0.5, lw=2, label='Deaths')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


