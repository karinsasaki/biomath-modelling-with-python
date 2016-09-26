# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            Dynamics of transcription and translation
on.
"""


# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline


# Parameters
kts1 = 0.0005
kex1 = 0.0037
krd1 = 8.085
ktl1 = 0.0205
kpd1 = 0.000125
#kts2 = 0.00045
#kex2 = 0.00005
#krd2 = 0.0937
#ktl2 = 0.00125
#kpd2 = 0.00014
params = (kts1, kex1, krd1, ktl1, kpd1)
#params = (kts2, kex2, krd2, ktl2, kpd2)

# Total numbers
Vratio = 0.5

# Initial conditions
mRNAnuc1 = 4*np.power(10,-3)
mRNAcyt1 = 1.06*np.power(10,-3)
Protein1 = 1.7*np.power(10,-6)
#mRNAnuc2 = 27*np.power(10,-2)
#mRNAcyt2 = 0.2*np.power(10,-3)
#Protein2 = 1.27*np.power(10,-3)
ini = (mRNAnuc1, mRNAcyt1, Protein1)
#ini = (mRNAnuc2, mRNAcyt2, Protein2)

# Time
dt = 0.001
t = np.arange(0,10,dt)

# ODE system (in format used for odeint)
def myODE(init,t,params):
    
    kts1, kex1, krd1, ktl1, kpd1 = params
    #kts2, kex2, krd2, ktl2, kpd2 = params
    
    mRNAnuc1, mRNAcyt1, Protein1 = init[0], init[1], init[2]
    #mRNAnuc2, mRNAcyt2, Protein2 = init[0], init[1], init[2]

    # rates
    vts = kts1*Hog1P2nuc
    vex = kex1*mRNAnuc1
    vrd = kts1*mRNAcyt1
    vtl = ktl1*mRNAcyt1
    vpd = kpd1*Protein1
    #vts = kts2*Hog1P2nuc
    #vex = kex2*mRNAnuc2
    #vrd = kts2*mRNAcyt2
    #vtl = ktl2*mRNAcyt2
    #vpd = kpd2*Protein2    
    
    # equations
    dmRNAnuc1 = vts - vex
    dmRNAcyt1 = vex*(Vnuc/Vcyt) - vrd - mRNAcyt1*Vratio
    dProtein1 = vtl - vpd - Protein1*Vratio
    #dmRNAnuc2 = vts - vex
    #dmRNAcyt2 = vex*(Vnuc/Vcyt) - vrd - mRNAcyt2*Vratio
    #dProtein2 = vtl - vpd - Protein2*Vratio    
    
    # conservation relations
    
    
    return (dmRNAnuc1, dmRNAcyt1, dProtein1)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 
tmRNAnuc1 = solution[:,0]
tmRNAcyt1 = solution[:,1]
tProtein1 = solution[:,2]
#tmRNAnuc1 = solution[:,0]
#tmRNAcyt1 = solution[:,1]
#tProtein1 = solution[:,2]


# Show over time
fig, ax = plt.subplots()
ax.plot(t, tmRNAnuc1, 'r', label='mRNAnuc1')
ax.plot(t, tmRNAcyt1, 'g', label='mRNAcyt1')
ax.plot(t, tProtein1, 'y', label='Protein1')
legend = ax.legend(loc='center', shadow=True)
plt.show()