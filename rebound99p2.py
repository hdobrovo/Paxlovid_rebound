import numpy as np
eps=0
CF=0.31
b,k,d,p,c =4.71*10**(-8)/CF,5.0,1.07,3.07*CF,2.4
T0,E0, I0, V0=4e8,0,0,1
b,k,d,p,c =4.71*10**(-8)/CF,5.0,1.07,3.07*CF,2.4
T1,E1, I1, V1=4e8,0,0,1
tau=0.01 
Max=30 
INPUT = np.array((T0,E0, I0, V0))
INPUT2 = np.array((T1,E1, I1, V1))
def stoc_eqs_tauleap(INP):
    X=INP
    Rate=np.zeros((10))
    Transitions=np.zeros((5,4))
    Rate[0] = b*X[0]*X[3]; Transitions[0,:]=([-1, +1, 0, 0]);
    Rate[1] = k*X[1];  Transitions[1,:]=([0, -1, +1, 0]);
    Rate[2] = d*X[2];  Transitions[2,:]=([0, 0, -1, 0]);
    Rate[3] = (1-eps)*p*X[2];  Transitions[3,:]=([0, 0, 0, +1]);
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
  

    for tt in time:
        result=stoc_eqs_tauleap(INPUT)
        T.append(result[0])
        E.append(result[1])
        I.append(result[2])
        V.append(result[3])
        INPUT=result

    return [T, E, I, V]

count=0
data=np.zeros((100,100))
n_simulations=100

for TT in np.arange(3,10,.07):
        count1=0
        for DELAY in np.arange(.5,7,.065):
            rebound=0
            for j in range (1,n_simulations+1):
           
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
                 rebound=rebound+1
            
           
            data[count1,count]=rebound
            count1=count1+1
        count=count+1
        
np.savetxt('rebound99p2.dat',data)        
        
