#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:27:20 2020

@author: leonardo
"""

import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from lmfit import minimize, Parameters, Parameter, report_fit, Model
from scipy import interpolate

np.random.seed(123)



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








N=83e6
I0 = 10
E0 = I0
A0 = I0
R0 = 0
D0 = 0
P0 = 0
Y0 = 0

S0=N-E0-A0-I0-R0-D0-P0


lb=0
ub=1
params = Parameters()

params = Parameters()
params.add('alpha0', value=np.random.uniform(1e-3,1e-1), min=lb, max=1)
params.add('alpha1', value=np.random.uniform(1e-1,1), min=lb, max=1)
params.add('beta', value= np.random.uniform(.8,2), min=lb, max=3)
params.add('gamma', value= np.random.uniform(1e-1,1), min=lb, max=1)
params.add('lambda0', value= np.random.uniform(0,0.03), min=0, max=1)
params.add('lambda1', value= np.random.uniform(0,0.03), min=0, max=1)
params.add('kappa0', value=np.random.uniform(0,0.03), min=0, max=1)
params.add('kappa1', value= np.random.uniform(0,0.03), min=0, max=1)
params.add('tau0', value= np.random.uniform(1e-3,1e-1), min=lb, max=1)
params.add('tau1', value= np.random.uniform(1e-1,1), min=lb, max=ub)
params.add('sigma', value= np.random.uniform(0.4,0.7), min=lb, max=1)

dt=1/24
tf = 180
tl = int(tf/dt)
t = np.linspace(0,tf,tl)
t2= np.linspace(t[0],t[len(t)-1],int((t[len(t)-1]-t[0])/10))
timesimul=pd.to_datetime(np.arange(pd.to_datetime('2021-01-01'), pd.to_datetime('2021-06-01'), dtype='datetime64[D]'))

S = np.zeros(tl)
E = np.zeros(tl)
A = np.zeros(tl)
I = np.zeros(tl)
R = np.zeros(tl)
D = np.zeros(tl)
P = np.zeros(tl)
Y = np.zeros(tl)
u = [S0,E0,A0,I0,R0,D0,P0,Y0]

#sir_out = pd.DataFrame(simulate(t,u,params)) #In case you want to run only 1 sumulation
sir_out = pd.DataFrame(simulate_N_times(t, u, params))

Yinter=Sys_interpolation(sir_out,t,t2)

plt.figure(1)
plt.style.use("ggplot")
#sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=2)
iline = plt.plot("t","E","",data=sir_out,color="cyan",linewidth=2)
iline = plt.plot("t","I","",data=sir_out,color="green",linewidth=2)
iline = plt.plot("t","A","",data=sir_out,color="gray",linewidth=2)
#rline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=2)
rline = plt.plot("t","D","",data=sir_out,color="black",linewidth=2)
plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")
plt.title('No interpolation')
legend = plt.legend(title="Population")#,loc=5,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)


plt.figure(2)
plt.style.use("ggplot")
#sline = plt.plot("t","S","",data=Yinter,color="red",linewidth=2)
iline = plt.plot("t","E","",data=Yinter,color="cyan",linewidth=2)
iline = plt.plot("t","I","",data=Yinter,color="green",linewidth=2)
iline = plt.plot("t","A","",data=Yinter,color="gray",linewidth=2)
#rline = plt.plot("t","R","",data=Yinter,color="blue",linewidth=2)
rline = plt.plot("t","D","",data=Yinter,color="black",linewidth=2)
plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")
plt.title('Interpolation')
legend = plt.legend(title="Population")#,loc=5,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)






