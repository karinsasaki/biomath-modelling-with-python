# -*- coding: utf-8 -*-

# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#%matplotlib inline

# Parameters
k0 = 1.2
k1 = 0.96
k2 = 1.18
k3 = 1
k4 = 1
k5 = 1
k6 = 1

params = (k0, k1, k2, k3, k4, k5, k6)

# Initial conditions
Receptor_0 = 1
Phosphotase_0 = 1
MAPKKKP_0 = 0
MAPKKP_0 = 0
MAPKP_0 = 0
MAPKKK_0 = 1
MAPKK_0 = 1
MAPK_0 = 1

#Receptor_0 = 1
#Phosphotase_0 = 1
#MAPKKK_0 = 5.306*np.power(10,-3)
#MAPKKKP_0 = 0.0067 - 5.306*np.power(10,-3)
#MAPKK_0 = 42.9*np.power(10,-3)
#MAPKKP_0 = 0.053 - (42.9*np.power(10,-3)+8.38*np.power(10,-3))
#MAPK_0 = 79*np.power(10,-3) 
#MAPKP_0 = 0.9*np.power(10,-3) + 0.03


ini = (Receptor_0, Phosphotase_0, MAPKKKP_0, MAPKKP_0, MAPKP_0, MAPKKK_0, MAPKK_0, MAPK_0)

# time
dt = 0.0001
t = np.arange(0,10,dt)


# ODE system (in format used for odeint)
def myODE(init,t,params):
    
    k0, k1, k2, k3, k4, k5, k6 = params
    
    Receptor, Phosphotase, MAPKKKP, MAPKKP, MAPKP, MAPKKK, MAPKK, MAPK = init

    # rates
    v0 = k0*Receptor
    v1 = k1*MAPKKK*Receptor
    v2 = k2*MAPKKKP*Phosphotase
    v3 = k3*MAPKK*MAPKKKP    
    v4 = k4*MAPKKP*Phosphotase
    v5 = k5*MAPK*MAPKKP    
    v6 = k6*MAPKP*Phosphotase

    # equations
    dReceptor = -v0
    dMAPKKKP = v1 - v2
    dMAPKKP = v3 - v4
    dMAPKP = v5 - v6
    dMAPKKK = -dMAPKKKP
    dMAPKK = -dMAPKKP
    dMAPK = -dMAPKP   
    dPhosphotase = 0
        
    return (dReceptor, dPhosphotase, dMAPKKKP, dMAPKKP, dMAPKP, dMAPKKK, dMAPKK, dMAPK)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 
tReceptor = solution[:,0]
tPhosphatase = solution[:,1]
tMAPKKKP = solution[:,2]
tMAPKKP = solution[:,3]
tMAPKP = solution[:,4]
tMAPKKK = solution[:,5]
tMAPKK = solution[:,6]
tMAPK = solution[:,7]

# Show over time
fig, ax = plt.subplots()
#ax.plot(t, tReceptor, label='Receptor')
ax.plot(t, tMAPKKKP, label='MAPKKKP')
ax.plot(t, tMAPKKP, label='MAPKKP')
ax.plot(t, tMAPKP, label='MAPKP')
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()