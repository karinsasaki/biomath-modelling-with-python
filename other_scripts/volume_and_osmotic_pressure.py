# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            Dynamics of the of the volume and osmotic pressure
on.
"""

# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
#%matplotlib inline


# Parameters
lambdam = 0.3
Lp = 1.19*math.pow(10,-12)
Vtotal_0 = 6.5*math.pow(10,-17)
Vb = 0.4*Vtotal_0
Vnuc = 0.15*Vtotal_0
Vcyt_0 = Vtotal - Vnuc
VPIt0 = 0.63*Vos0
G = 7.85*math.pow(10,-11)
PIt0 = 0.875*math.pow(10,6)
PIe0 = 0.625*math.pow(10,6)

params = (lambdam, Lp, Vb, Vnuc, Vcyt_0, VPIt0, VOS0, G, c0, PIt0, PIi0, PIe0)

# Total numbers

# Initial conditions
PIi0 = 1.5*math.pow(10,6)
c0 = 600
Vos0 = Vtotal0 - Vb

ini = (PIi0, c0, Vos0)

# Time
dt = 0.001
t = np.arange(0,10,dt)


# ODE system (in format used for odeint)
def myODE(init,t,params):
    
     lambdam, Lp, Vb, Vnuc, Vcyt_0, VPIt0, VOS0, G, c0, PIt0, PIi0, PIe0 = params
    
     PIi, c, Vos = init[0], init[1], init[2]
    
    
    # rates
    PIe = PIfe*(1 - math.exp(-(t-t_0)/lambdam) + PIe_0*math.exp(-(t-t_0)/lambdam))
    PI_t = PIt_0*(1-(Vos_0-Vos)/Vos_0-VPI_t0) ##KS
    
    # equations
    dPIi = PIi*dc/c
    dc = v12 - v13 - c*Vratio
    dVos = -G*Lp*(PIt+PIe-PIi)
    
    # conservation relations
    Vos = Vtotal - Vb
    dVtotal = dVcyt #(=dVos)
    Vratio = dVos/Vos
    
    return (dPIi, dc, dVos)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 


# Show over time
fig, ax = plt.subplots()


legend = ax.legend(loc='center', shadow=True)
plt.show()
