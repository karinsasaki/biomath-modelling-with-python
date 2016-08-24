# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            Dynamics of the phosphorelay module
on.
"""


# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline


# Parameters
k1TCS = 5
k2TCS = 50
k2mTCS = 50
k3TCS = 50
k4TCS = 0.415
params = (k1TCS, k2TCS, k2mTCS, k3TCS, k4TCS)

# Total numbers
Sln1T = 0.016
Ypd1T = 0.156
Ssk1T = 0.029
V = 1 ##KS

# Initial conditions
Sln1P = 2.25*np.power(10,-3)
Ypd1P = 36*np.power(10,-3)
Ssk1P = 1.88*np.power(10,-3)
Sln1 = Sln1T - Sln1P ##KS
Ypd1 = Ypd1T - Ypd1P ##KS ##KS
Ssk1 = Ssk1T - Ssk1P ##KS
ini = (Sln1, Sln1P, Ypd1, Ypd1P, Ssk1, Ssk1P)

# Time
dt = 0.001
t = np.arange(0,10,dt)

# ODE system (in format used for odeint)
def myODE(init,t,params):
    
    k1TCSTCS, k2TCS, k2mTCS, k3TCS, k4TCS = params
    
    Sln1, Sln1P, Ypd1, Ypd1P, Ssk1, Ssk1P = init[0], init[1], init[2], init[3], init[4], init[5]

    # rates
    v1TCS = k1TCS*Sln1
    v2TCS = k2TCS*Sln1P*Ypd1 - k2mTCS*Sln1*Ypd1P
    v3TCS = k3TCS*Ssk1*Ypd1P
    v4TCS = k4TCS*Ssk1P
    
    # equations
    dSln1 = -v1TCS + v2TCS - Sln1*V
    dSln1P = v1TCS - v2TCS - Sln1P*V
    dYpd1 = -v2TCS + v2TCS - Ypd1*V
    dYpd1P = v2TCS - v3TCS - Ypd1P*V
    dSsk1 = -v3TCS + v4TCS - Ssk1*V
    dSsk1P = v3TCS - v4TCS - Ssk1P*V
    
    # conservation relations
    Sln1 = Sln1T - Sln1P ##KS
    Ypd1 = Ypd1T - Ypd1P ##KS ##KS
    Ssk1 = Ssk1T - Ssk1P ##KS
    
    return (dSln1, dSln1P, dYpd1, dYpd1P, dSsk1, dSsk1P)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 
tSln1 = solution[:,0]
tSln1P = solution[:,1]
tYpd1 = solution[:,2]
tYpd1P = solution[:,3]
tSsk1 = solution[:,4]
tSsk1P =solution[:,5]

# Show over time
fig, ax = plt.subplots()
ax.plot(t, tSln1, 'r', label='Sln1')
ax.plot(t, tSln1P, 'g', label='Sln1P')
ax.plot(t, tYpd1, 'y', label='Ypd1')
ax.plot(t, tYpd1P, 'c', label='Ypd1P')
ax.plot(t, tSsk1, 'k', label='Ssk1')
ax.plot(t, tSsk1P, 'b', label='Ssk1P')
legend = ax.legend(loc='center', shadow=True)
plt.show()