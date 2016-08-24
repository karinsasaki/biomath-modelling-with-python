# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            Dynamics of the MAP kinase cascade
"""




# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib inline


# Parameters
k1MAP = 1.438
k2MAP = 1.438
k3MAP = 1.438
k4MAP = 1.438
k5MAP = 1.438
k1mMAP = 0.011
k2mMAP = 0.011
k2mMAP = 0.011
k3mMAP = 0.011
k4mMAP = 0.011
k5mMAP = 0.011
kHog1P2trans = 0.029
kHog1dephos = 0.0053
kHog1trans1 = 0.110
kHog1trans2 = 0.091
params = (k1MAP, k2MAP, k3MAP, k4MAP, k1mMAP, k2mMAP, k3mMAP, k4mMAP, kHog1P2trans, kHog1dephos, kHog1trans1, kHog1trans2)

# Total numbers
Ssk2T = 0.0067
Pbs2T = 0.053
Hog1T = 0.167

Vcyt = 1 ##KS
Vnuc = 0.5 ##KS
Vratio = 0.5 ##KS
VT = 1.5 ##KS
protein2 = 1 ##KS
protein2_0=1 ##KS

# Initial conditions
Ssk2 = 5.306*np.power(10,-3)
Pbs2 = 42.9*np.power(10,-3)
Pbs2P = 8.38*np.power(10,-3)
Hog1 = 79*np.power(10,-3)
Hog1P = 15.5*np.power(10,-3)
Hog1P2cyt = 0.9*np.power(10,-3)
Hog1P2nuc = 0.03
Hog1nuc = 0.5

Ssk2P = Ssk2T - Ssk2 ##KS
Pbs2P2 = Pbs2T - (Pbs2 + Pbs2P) ##KS
Hog1T = (Hog1 + Hog1P + Hog1P2cyt)*(Vcyt/VT) + (Hog1P2nuc + Hog1nuc)*(Vnuc/VT) ##KS

ini = (Ssk2, Ssk2P, Pbs2, Pbs2P, Pbs2P2, Hog1, Hog1P, Hog1P2cyt, Hog1P2nuc, Hog1nuc)

# Time
dt = 0.001
t = np.arange(0,10,dt)

# ODE system (in format used for odeint)
def myODE(init,t,params):
    
    k1MAP, k2MAP, k3MAP, k4MAP, k1mMAP, k2mMAP, k3mMAP, k4mMAP, kHog1P2trans, k1Hog1dephos, k1Hog1trans1, k1Hog1trans2 = params
    
    Ssk2, Ssk2P, Pbs2, Pbs2P, Pbs2P2, Hog1, Hog1P, Hog1P2cyt, Hog1P2nuc, Hog1nuc = init[0], init[1], init[2], init[3], init[4], init[5], init[6], init[7], init[8], init[9]

    # rates
    v1MAP = k1MAP*Ssk2*1#Ssk1 ##KS
    v1mMAP = k1mMAP*Ssk2P
    v2MAP = k2MAP*Pbs2*Ssk2P
    v2mMAP = k2mMAP*Pbs2P
    v3MAP = k3MAP*Pbs2P*Ssk2P
    v3mMAP = k3mMAP*Pbs2P2
    v4MAP = k4MAP*Hog1*Pbs2P2
    v4mMAP = k4mMAP*Hog1P
    v5MAP = k5MAP*Hog1P*Pbs2P2
    v5mMAP = k5mMAP*Hog1P2cyt
    vHog1trans1 = kHog1trans1*Hog1nuc
    vHog1trans2 = kHog1trans2*Hog1
    vHog1P2trans = kHog1P2trans*Hog1P2nuc
    vHog1dephos = (protein2/protein2_0)*kHog1dephos*Hog1P2nuc
    Hog1P2 = Hog1P2cyt*(Vcyt/VT)+Hog1P2nuc*(Vnuc/VT)
    
    # equations
    dSsk2 = -v1MAP + v1mMAP - Ssk2*Vratio
    dSsk2P = v1MAP - v1mMAP - Ssk2P*Vratio
    dPbs2 = -v2MAP + v2mMAP - Pbs2*Vratio
    dPbs2P = v2MAP - v2mMAP - v2MAP + v3mMAP - Pbs2P*Vratio
    dPbs2P2 = v3MAP - v3mMAP - Pbs2P2*Vratio
    dHog1 = -v4MAP + v4mMAP - vHog1trans2 + vHog1trans1*(Vnuc/Vcyt) - Hog1*Vratio
    dHog1P = v4MAP - v4mMAP - v5MAP + v5mMAP - Hog1P*Vratio
    dHog1P2cyt = v5MAP - v5mMAP - vHog1P2trans - Hog1P2cyt*Vratio
    dHog1P2nuc = vHog1P2trans*(Vcyt/Vnuc) - vHog1dephos
    dHog1nuc = vHog1trans2*(Vcyt/Vnuc) - vHog1trans1 + vHog1dephos
    
    # conservation relations
    Ssk2T = Ssk2 + Ssk2P ##KS 
    Pbs2T = Pbs2 + Pbs2P + Pbs2P2 ##KS
    Hog1T = (Hog1 + Hog1P + Hog1P2cyt)*(Vcyt/VT) + (Hog1P2nuc + Hog1nuc)*(Vnuc/VT) ##KS
    
    return (dSsk2, dSsk2P, dPbs2, dPbs2P, dPbs2P2, dHog1, dHog1P, dHog1P2cyt, dHog1P2nuc, dHog1nuc)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 
tSsk2 = solution[:,0]
tSsk2P = solution[:,1]
tPbs2 = solution[:,2]
tPbs2P = solution[:,3]
tPbs2P2 = solution[:,4]
tHog1 = solution[:,5]
tHog1P = solution[:,6]
tHog1P2cyt = solution[:,7]
tHog1P2nuc = solution[:,8]
tHog1nuc = solution[:,9]

# Show over time
fig, ax = plt.subplots()
ax.plot(t, tSsk2, label='Ssk2')
ax.plot(t, tSsk2P, label='Ssk2P')
ax.plot(t, tPbs2, label='Pbs2')
ax.plot(t, tPbs2P, label='Pbs2P')
ax.plot(t, tPbs2P2, label='Pbs2P2')
ax.plot(t, tHog1, label='Hog1')
ax.plot(t, tHog1P, label='Hog1P')
ax.plot(t, tHog1P2cyt, label='Hog1P2cyt')
ax.plot(t, tHog1P2nuc, label='Hog1P2nuc')
ax.plot(t, tHog1nuc, label='Hog1nuc')

legend = ax.legend(loc='center', shadow=True)
plt.show()