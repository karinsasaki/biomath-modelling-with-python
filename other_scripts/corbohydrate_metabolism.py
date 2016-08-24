# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            Dynamics of the of the carbohydrate metabolism
on.
"""

# import modules and libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline


# Parameters
v1_0 = 1.296
k2 = 1.777
k21 = 0.2
k22 = 1.2
k23 = 0.2
k3 = 0.895
k31 = 0.01
k32 = 0.012
k4 = 1.8764*np.power(10,-3)
Keq4 = 0.81
kFBP4 = 0.054
kDHAP4 = 2
kGAP4 = 2
vbf4 = 5
KIGAP4 = 10
k5 =  190
Keq5 = 0.45
KDHAP5 = 0.38
KGAP5 = 0.064
k6 = 45.127
k7 = 639.137
KPyr7 = 70
KNAD7 = 160
KI7 = 20
k8 = 5.6425
k81  = 1.2
k82 = 0.6
k9 = 0.8090
KPyr9 = 0.92
KATP9 = 13.2
k10 = 1.9377
KG6P10 = 2.7
KATP10 = 0.4
k11 = 7.1507
kd11 = 0.0037
kn11 = 0.6
kh11 = 0.2
Keq11 = 37
k12 = 0.0162
k13 = 0.005 ##KS
kx = 6.67*np.power(10,-10)
n13 = 4
k14 = 384.024
KADP14 = 0.42
k15 = 1.1235
k16 = 99.887
KATP16 = 5.0

params = (v1_0, k2, k21, k22, k23, k3, k31, k32, k4, Keq4, kFBP4, kDHAP4, kGAP4, vbf4, KIGAP4, k5, Keq5, KDHAP5, KGAP5, k6, k7, KPyr7, KNAD7, KI7, k8, k81, k82, k9, KPyr9, KATP9, k10, KG6P10, KATP10, k11, kd11, kn11, kh11, Keq11, k12, k13_0, kx, n13, k14, KADP14, k15, k16, KATP16)

# Total numbers
Vratio = 0.5

# Initial conditions
Gluc = 1
G6P = 0.3
FBP = 0.11
GAP = 0.063
DHAP = 1.41
Pyr = 0.79
G3P = 0.063
Glycin = 0.2
Glycex = 0.2
NADH = 0.04
NADp = 1.03 ##KS
ATP = 2.1
ADP = 0.6
Protein1_0 = 1 ##KS

ini = (Gluc, G6P, FBP, GAP, DHAP, Pyr, G3P, Glycin, Glycex, NADH, NADp, ATP, ADP, Protein1_0)

# Time
dt = 0.001
t = np.arange(0,10,dt)

# ODE system (in format used for odeint)
def myODE(init,t,params):
    v1_0, k2, k21, k22, k23, k3, k31, k32, k4, Keq4, kFBP4, kDHAP4, kGAP4, vbf4, KIGAP4, k5, Keq5, KDHAP5, KGAP5, k6, k7, KPyr7, KNAD7, KI7, k8, k81, k82, k9, KPyr9, KATP9, k10, KG6P10, KATP10, k11, kd11, kn11, kh11, Keq11, k12, k13, kx, n13, k14, KADP14, k15, k16, KATP16 = params
    
    Gluc, G6P, FBP, GAP, DHAP, Pyr, G3P, Glycin, Glycex, NADH, NADp, ATP, ADP, Protein1_0 = init[0], init[1], init[2], init[3], init[4], init[5], init[6], init[7], init[8], init[9], init[10], init[11], init[12], init[13]
    
    Glyc = Glycin + Glycex ##KS
    
    # rates
    print 'rates'
    v1 = (Protein1/Protein1_0)*v1_0
    v2 = (Protein1/Protein1_0)*(k2/(1+(k21/ATP)*(1+(ADP/k22))+(k23/Gluc)+(k21/ATP)*(k23/Gluc)*(1+(ADP/k22))))
    v3 = k3*(ATP/(k31 + ATP))*(G6P/(k32 + G6P))
    v4 = k4*((FBP-(GAP*DHAP/Keq4)))/(kFBP4 + FBP + ((kDHAP4*GAP)/(Keq4*vbf4)) + kGAP4/(Keq4*vbf4) + FBP*GAP/KIGAP4 + (GAP*DHAP)/(Keq4*vbf4))
    v5 = k5*((DHAP - (GAP/Keq5))/(KDHAP5*(1+(GAP/KGAP5)+DHAP)))
    v6 = k6*GAP*NAD*ADP
    v7 = k7*((ADP*Pyr*NAD)/(KNAD7*Pyr + KPyr7*NAD + KI7*NADH))
    v8 = k8*(Pyr/(k81+Pyr))*(NAD/(k82+NAD))
    v9 = k9*(Pyr/(KPyr9+Pyr))*(ATP/(KATP9+ATP))
    v10 = k10*(G6P/(KG6P10+G6P))*(ATP*(KATP10+ATP))
    v11 = (Protein1/Protein1_0)*k11*((DHAP*NADH-((G3P*NAD)/Keq11))/(1+(DHAP/kd11)+(NADH/kn11)+((G3P*NAD)/kh11)))
    v12 = (Protein1/Protein1_0)*k12*G3P
    v13 = k13*(Glyc-Glycex)
    v14 = k14*(NADH*ADP/(KADP14+ADP))
    v15 = k15*NAD
    v16 = k16*(ATP/(KATP16+ATP))

    # equations
    print 'equations'
    dGluc = v1 - v2 - Gluc*Vratio
    dG6P = v2 - v3 - v10 - G6P*Vratio
    dFBP = v3 - v4 - FBP*Vratio
    dGAP = v4 + v5 - v6 - GAP*Vratio
    dDHAP = v4 - v5 - v11 - DHAP*Vratio
    dPyr = v6 - v7 - v8 - v9 - Pyr*Vratio
    dG3P = v11 - v12 - G3P*Vratio
    dGlycin = v12 - v13 - Glyc*Vratio
    dGlycex = v13/kx
    dNADH = v6 + 4*v7 + v8 + v11 + v14 - v15 - NADH*Vratio
    dNADp = -v6 - 4*v7 - v8 - v11 - v14 + v15 - NADp*Vratio
    dATP = -v2 - v3 + 2*v6 + v7 - v9 - v10 + 3*v14 - v16 - ATP*Vratio
    dADP = v2 + v3 - 2*v6 - v7 + v9 + v10 - 3*v14 + v16 - ADP*Vratio

    # conservation relations

    return (dGluc, dG6P, dFBP, dGAP, dDHAP, dPyr, dG3P, dGlycin, dGlycex, dNADH, dNADp, dATP, dADP, Protein1_0)

# Solve using odeint
solution = odeint(myODE,ini,t,args=(params,)) 
tGluc = solution[:,0]
tG6P = solution[:,1]
tFBP = solution[:,2]
tGAP = solution[:,3]
tDHAP = solution[:,4]
tPyr = solution[:,5]
tG3P = solution[:,6]
tGlycin = solution[:,7]
tGlycex = solution[:,8]
tNADH = solution[:,9]
tNAD = solution[:,10]
tATP = solution[:,11]
tADP = solution[:,12]
tProtein1_0 = solution[:,13] ##KS

# Show over time
fig, ax = plt.subplots()
ax.plot(t, tGluc, label='Gluc')
ax.plot(t, tG6P, label='G6P')
ax.plot(t, tFBP, label='FBP')
ax.plot(t, tGAP, label='GAP')
ax.plot(t, tDHAP, label='DHAP')
ax.plot(t, tPyr, label='Pyr')
ax.plot(t, tG3P, label='G3P')
ax.plot(t, tGlycin, label='Glycin')
ax.plot(t, tGlycex, label='Glycex')
ax.plot(t, tNADH, label='NADH')
ax.plot(t, tNAD, label='NAD')
ax.plot(t, tATP, label='ATP')
ax.plot(t, tADP, label='ADP')
ax.plot(t, tProtein1_0, label='Protein1_0') ##KS

legend = ax.legend(loc='center', shadow=True)
plt.show()
