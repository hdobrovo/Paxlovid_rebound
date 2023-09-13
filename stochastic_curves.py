import numpy as np
import matplotlib.pyplot as plt
import random
import pylab as pl
import seaborn as sns; sns.set(color_codes=True)
import datetime
import pandas as pd
eps=0
CF=0.31
b,k,d,p,c =4.71*10**(-8)/CF,5.0,1.07,3.07*CF,2.4
T0,E0, I0, V0=4e8,0,0,1
b,k,d,p,c =4.71*10**(-8)/CF,5.0,1.07,3.07*CF,2.4
T1,E1, I1, V1=4e8,0,0,1
tau=0.01 #dt1/0.31
Max=30 #maximum days of infection
INPUT = np.array((T0,E0, I0, V0))
INPUT2 = np.array((T1,E1, I1, V1))
def stoc_eqs_tauleap(INP):
    X=INP
    Rate=np.zeros((10))##propensity function [b1*T*V1 k1*E1 d1*I1 p1*I1 c1*V1 b2*T*V2 k2*E2 d2*I2 p2*I2 c2*V2]
    Transitions=np.zeros((5,4))##stoichiometric matrix, each row of which is a transition vector
    Rate[0] = (1-eps)*b*X[0]*X[3]; Transitions[0,:]=([-1, +1, 0, 0]);
    Rate[1] = k*X[1];  Transitions[1,:]=([0, -1, +1, 0]);
    Rate[2] = d*X[2];  Transitions[2,:]=([0, 0, -1, 0]);
    Rate[3] = p*X[2];  Transitions[3,:]=([0, 0, 0, +1]);
    Rate[4] = c*X[3];  Transitions[4,:]=([0, 0, 0, -1]);
    for i in range(5):
        leap=np.random.poisson(Rate[i]*tau);

        if np.size(np.where(Transitions[i,:]<0)[0])==0:
           Use=leap
        else:
           Use=min([leap, X[np.where(Transitions[i,:]<0)[0][0]]]);
        X=X+Transitions[i,:]*Use;
    return X

def Stoch_Iteration(INPUT,time):

    T=[INPUT[0]]
    E=[INPUT[1]]
    I=[INPUT[2]]
    V=[INPUT[3]]
    #np.array((T,E,I,V))

    for tt in time:
        result=stoc_eqs_tauleap(INPUT)# create second graph
        T.append(result[0])
        E.append(result[1])
        I.append(result[2])
        V.append(result[3])
        INPUT=result#adding in to teh stoch info thats already there

    return [T, E, I, V]



#plt.figure(figsize=(12, 9))
count=0
data=np.zeros((100,100))
n_simulations=10


TT=10
DELAY=6
rebound = 0
plt.figure(1)
plt.figure(2)
colors=['b','g','r','c','m','k','tab:orange','tab:purple','tab:brown','tab:pink']
for j in range (1,n_simulations+1):
     print('simulation=',j)
     eps=0

     time=np.arange(0,DELAY,.01)
     time2=np.arange(DELAY,DELAY+TT,.01)
     time3=np.arange(DELAY+TT,30,.01)
     [T,E,I,V]=Stoch_Iteration(INPUT,time)
     eps=0.99
     INPUT2=[T[-1],E[-1],I[-1],V[-1]]
     [T1,E1,I1,V1]=Stoch_Iteration(INPUT2,time2)
     eps=0
     INPUT3=[T1[-1],E1[-1],I1[-1],V1[-1]]
     [T2,E2,I2,V2]=Stoch_Iteration(INPUT3,time3)
     t=np.array(time)
     r=V1[-1]+V1[-1]/10
     if max(V2)>V1[-1]+V1[-1]/10:
         print('rebound')
         rebound=rebound+1
     plt.figure(1)    
     plt.semilogy(time,V[:-1],colors[j-1],linewidth=2)
     plt.semilogy(time2,V1[:-1],colors[j-1],linewidth=2)
     plt.semilogy(time3,V2[:-1],colors[j-1],linewidth=2)
     plt.figure(2)
     plt.semilogy(time,I[:-1],colors[j-1],linewidth=2)
     plt.semilogy(time2,I1[:-1],colors[j-1],linewidth=2)
     plt.semilogy(time3,I2[:-1],colors[j-1],linewidth=2)     

count=count+1
plt.figure(1)
plt.ylabel('Viral titer [V]', fontsize=30)
plt.xlabel('Time (d)', fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=25)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.xlim(0,30)
plt.savefig('virus_b_10_6.pdf', bbox_inches='tight', pad_inches=0.2)
plt.figure(2)
plt.ylabel('Infected cells', fontsize=30)
plt.xlabel('Time (d)', fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=25)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.xlim(0,30)
plt.savefig('cells_b_10_6.pdf', bbox_inches='tight', pad_inches=0.2)
