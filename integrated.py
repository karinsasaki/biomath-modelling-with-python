# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 2016

@author:    Karin Sasaki

@descript:  Mathematical model of the High Osmolarity Gycerol (HOG) system of S. cerevisiae
            All modules integrated
on.
"""

# ------------------------------------
# import modules and libraries
# ------------------------------------
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from mpmath import *
#%matplotlib inline

# ------------------------------------
# Parameters
# ------------------------------------

    # phosphorelay module
k1TCS0 = 5
n1 = 2
k2TCS = 50
k2mTCS = 50
k3TCS = 50
k4TCS = 0.415

    # MAP kinase cascade
k1MAP = 1.438
k2MAP = 1.438
k3MAP = 1.438
k4MAP = 1.438
k5MAP = 1.438
k1mMAP = 0.011
k2mMAP = 0.011
k3mMAP = 0.011
k4mMAP = 0.011
k5mMAP = 0.011
kHog1P2trans = 0.029
kHog1dephos = 0.0053
kHog1trans1 = 0.110
kHog1trans2 = 0.091

    # transcritpion and translation
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

    # carbohydrate metabolism
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
k13_0 = 0.005 ##KS
kx = 6.67*np.power(10,-10)
n13 = 4
k14 = 384.024
KADP14 = 0.42
k15 = 1.1235
k16 = 99.887
KATP16 = 5.0  

# volume and osmotic pressure  
lambdam = 0.3
Lp = 1.19*math.pow(10,-12)
G = 7.85*math.pow(10,-11)

params = (k1TCS0, n1, k2TCS, k2mTCS, k3TCS, k4TCS, k1MAP, k2MAP, k3MAP, k4MAP, k5MAP, k1mMAP, k2mMAP, k3mMAP, k4mMAP, k5mMAP, kHog1P2trans, kHog1dephos, kHog1trans1, kHog1trans2, kts1, kex1, krd1, ktl1, kpd1, v1_0, k2, k21, k22, k23, k3, k31, k32, k4, Keq4, kFBP4, kDHAP4, kGAP4, vbf4, KIGAP4, k5, Keq5, KDHAP5, KGAP5, k6, k7, KPyr7, KNAD7, KI7, k8, k81, k82, k9, KPyr9, KATP9, k10, KG6P10, KATP10, k11, kd11, kn11, kh11, Keq11, k12, k13_0, kx, n13, k14, KADP14, k15, k16, KATP16, lambdam, Lp, G)



# ------------------------------------
# Initial conditions
# ------------------------------------

    # volume and osmotic pressure
c_0 = 600

PIt_0 = 0.875*math.pow(10,6)
PIi_0 = 1.5*math.pow(10,6)
PIe_0 = 0.625*math.pow(10,6)

Vtotal_0 = 6.5*math.pow(10,-17)
Vb = 0.4*Vtotal_0
Vnuc = 0.15*Vtotal_0
Vcyt = Vtotal_0 - Vnuc
Vratio = Vcyt/Vnuc
Vos = Vtotal_0 - Vb
VPIt_0 = 0.63*Vos

#PIt_0 = PIi_0 - PIe_0


    # phosphorelay module
Sln1P_0 = 2.25*np.power(10,-3)
Ypd1P_0 = 36*np.power(10,-3)
Ssk1P_0 = 1.88*np.power(10,-3)
Sln1total = 0.016
Ypd1total = 0.156
Ssk1total = 0.029
Sln1_0 = Sln1total - Sln1P_0 ##KS
Ypd1_0 = Ypd1total - Ypd1P_0 ##KS ##KS
Ssk1_0 = Ssk1total - Ssk1P_0 ##KS


    # MAP kinase cascade
Ssk2_0 = 5.306*np.power(10,-3)
Pbs2_0 = 42.9*np.power(10,-3)
Pbs2P_0 = 8.38*np.power(10,-3)
Hog1_0 = 79*np.power(10,-3)
Hog1P_0 = 15.5*np.power(10,-3)
Hog1P2cyt_0 = 0.9*np.power(10,-3)
Hog1P2nuc_0 = 0.03
Hog1nuc_0 = 0.5

Ssk2total = 0.0067
Pbs2total = 0.053
Hog1total = 0.167

Ssk2P_0 = Ssk2total - Ssk2_0
Pbs2P2_0 = Pbs2total - (Pbs2_0 + Pbs2P_0)
Hog1T_0 = (Hog1_0 + Hog1P_0 + Hog1P2cyt_0)*(Vcyt/Vtotal_0) + (Hog1P2nuc_0 + Hog1nuc_0)*(Vnuc/Vtotal_0)


    # transcritpion and translation
mRNAnuc1_0 = 4*np.power(10,-3)
mRNAcyt1_0 = 1.06*np.power(10,-3)
Protein1_0 = 1.7*np.power(10,-6)
#mRNAnuc2_0 = 27*np.power(10,-2)
#mRNAcyt2_0 = 0.2*np.power(10,-3)
#Protein2_0 = 1.27*np.power(10,-3)


    # carbohydrate metabolism
Gluc_0 = 1
G6P_0 = 0.3
FBP_0 = 0.11
GAP_0 = 0.063
DHAP_0 = 1.41
Pyr_0 = 0.79
G3P_0 = 0.063
Glycin_0 = 0.2
Glycex_0 = 0.2
NADH_0 = 0.04
NADp_0 = 1.03 ##KS
ATP_0 = 2.1
ADP_0 = 0.6



ini = (c_0, PIt_0, PIi_0, PIe_0, Vtotal_0, Vb, Vnuc, Vcyt, Vratio, Vos, VPIt_0, Sln1P_0, Ypd1P_0, Ssk1P_0, Sln1total, Ypd1total, Ssk1total, Sln1_0, Ypd1_0, Ssk1_0, Ssk2_0, Pbs2_0, Pbs2P_0, Hog1_0, Hog1P_0, Hog1P2cyt_0, Hog1P2nuc_0, Hog1nuc_0, Ssk2total, Pbs2total, Hog1total, Ssk2P_0, Pbs2P2_0, Hog1T_0, mRNAnuc1_0, mRNAcyt1_0, Protein1_0, Gluc_0, G6P_0, FBP_0, GAP_0, DHAP_0, Pyr_0, G3P_0, Glycin_0, Glycex_0, NADH_0, NADp_0, ATP_0, ADP_0)


             
# ------------------------------------
# Time
# ------------------------------------

dt = 0.001
t = np.arange(0,10,dt)

# ------------------------------------
# ODE system (in format used for odeint)
# ------------------------------------

def myODE(init,t,params):
    k1TCS0, n1, k2TCS, k2mTCS, k3TCS, k4TCS, k1MAP, k2MAP, k3MAP, k4MAP, k5MAP, k1mMAP, k2mMAP, k3mMAP, k4mMAP, k5mMAP, kHog1P2trans, kHog1dephos, kHog1trans1, kHog1trans2, kts1, kex1, krd1, ktl1, kpd1, v1_0, k2, k21, k22, k23, k3, k31, k32, k4, Keq4, kFBP4, kDHAP4, kGAP4, vbf4, KIGAP4, k5, Keq5, KDHAP5, KGAP5, k6, k7, KPyr7, KNAD7, KI7, k8, k81, k82, k9, KPyr9, KATP9, k10, KG6P10, KATP10, k11, kd11, kn11, kh11, Keq11, k12, k13, kx, n13, k14, KADP14, k15, k16, KATP16, lambdam, Lp, G = params

    #c, PIt, PIi, PIe, Vtotal, Vb, Vnuc, Vcyt, Vratio, Vos, VPIt, Sln1P, Ypd1P, Ssk1P, Sln1total, Ypd1total, Ssk1total, Sln1, Ypd1, Ssk1, Ssk2, Pbs2, Pbs2P, Hog1, Hog1P, Hog1P2cyt, Hog1P2nuc, Hog1nuc, Ssk2total, Pbs2total, Hog1total, Ssk2P, Pbs2P2, Hog1T, mRNAnuc1, mRNAcyt1, Protein1, Gluc, G6P, FBP, GAP, DHAP, Pyr, G3P, Glycin, Glycex, NADH, NADp, ATP, ADP = init[0], init[1], init[2], init[3], init[4], init[5], init[6], init[7], init[8], init[9], init[10], init[11], init[12], init[13], init[14], init[15], init[16], init[17], init[18], init[19], init[20], init[21], init[22], init[23], init[24], init[25], init[26], init[27], init[28], init[29], init[30], init[31], init[32], init[33], init[34], init[35], init[36], init[37], init[38], init[39], init[40], init[41], init[42], init[43], init[44], init[45], init[46], init[47], init[48], init[49], init[50]
    c, PIt, PIi, PIe, Vtotal, Vb, Vnuc, Vcyt, Vratio, Vos, VPIt, Sln1P, Ypd1P, Ssk1P, Sln1total, Ypd1total, Ssk1total, Sln1, Ypd1, Ssk1, Ssk2, Pbs2, Pbs2P, Hog1, Hog1P, Hog1P2cyt, Hog1P2nuc, Hog1nuc, Ssk2total, Pbs2total, Hog1total, Ssk2P, Pbs2P2, Hog1T, mRNAnuc1, mRNAcyt1, Protein1, Gluc, G6P, FBP, GAP, DHAP, Pyr, G3P, Glycin, Glycex, NADH, NADp, ATP, ADP = init


    # **********
    # parameters
    # **********

    k1TCS = k1TCS0 * math.pow((PIt / PIt_0), n1)
    Vos_0 = 6.5 * math.pow(10, -17) - 0.4 * Vtotal_0
    Protein1_0 = 1.7 * np.power(10, -6)

    # **********
    # rates
    # **********
    
        # phosphorelay module
    v1TCS = k1TCS*Sln1
    v2TCS = k2TCS*Sln1P*Ypd1 - k2mTCS*Sln1*Ypd1P
    v3TCS = k3TCS*Ssk1*Ypd1P
    v4TCS = k4TCS*Ssk1P
    
        # MAP kinase cascade
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
    vHog1dephos = (Protein1/Protein1_0)*kHog1dephos*Hog1P2nuc
    Hog1P2 = Hog1P2cyt*(Vcyt/Vtotal)+Hog1P2nuc*(Vnuc/Vtotal)
            
        # transcritpion and translation
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
    
        # carbohydrate metabolism
    v1 = (Protein1/Protein1_0)*v1_0
    v2 = (Protein1/Protein1_0)*(k2/(1+(k21/ATP)*(1+(ADP/k22))+(k23/Gluc)+(k21/ATP)*(k23/Gluc)*(1+(ADP/k22))))
    v3 = k3*(ATP/(k31 + ATP))*(G6P/(k32 + G6P))
    v4 = k4*((FBP-(GAP*DHAP/Keq4)))/(kFBP4 + FBP + ((kDHAP4*GAP)/(Keq4*vbf4)) + kGAP4/(Keq4*vbf4) + FBP*GAP/KIGAP4 + (GAP*DHAP)/(Keq4*vbf4))
    v5 = k5*((DHAP - (GAP/Keq5))/(KDHAP5*(1+(GAP/KGAP5)+DHAP)))
    v6 = k6*GAP*NADp*ADP
    v7 = k7*((ADP*Pyr*NADp)/(KNAD7*Pyr + KPyr7*NADp + KI7*NADH))
    v8 = k8*(Pyr/(k81+Pyr))*(NADp/(k82+NADp))
    v9 = k9*(Pyr/(KPyr9+Pyr))*(ATP/(KATP9+ATP))
    v10 = k10*(G6P/(KG6P10+G6P))*(ATP*(KATP10+ATP))
    v11 = (Protein1/Protein1_0)*k11*((DHAP*NADH-((G3P*NADp)/Keq11))/(1+(DHAP/kd11)+(NADH/kn11)+((G3P*NADp)/kh11)))
    v12 = (Protein1/Protein1_0)*k12*G3P
    v13 = k13*(Glycin - Glycex)
    v14 = k14*(NADH*ADP/(KADP14+ADP))
    v15 = k15*NADp
    v16 = k16*(ATP/(KATP16+ATP))
    
        # volume and osmotic pressure
    
    # **********    
    # equations
    # **********
    
        # phosphorelay module
    dSln1 = -v1TCS + v2TCS - Sln1*Vratio    ##KSC
    dSln1P = v1TCS - v2TCS - Sln1P*Vratio   ##KSC
    dYpd1 = -v2TCS + v2TCS - Ypd1*Vratio   ##KSC
    dYpd1P = v2TCS - v3TCS - Ypd1P*Vratio   ##KSC
    dSsk1 = -v3TCS + v4TCS - Ssk1*Vratio   ##KSC
    dSsk1P = v3TCS - v4TCS - Ssk1P*Vratio   ##KSC
    
        # MAP kinase cascade
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
            
        # transcritpion and translation
    dmRNAnuc1 = vts - vex
    dmRNAcyt1 = vex*(Vnuc/Vcyt) - vrd - mRNAcyt1*Vratio
    dProtein1 = vtl - vpd - Protein1*Vratio
    #dmRNAnuc2 = vts - vex
    #dmRNAcyt2 = vex*(Vnuc/Vcyt) - vrd - mRNAcyt2*Vratio
    #dProtein2 = vtl - vpd - Protein2*Vratio
        
        # carbohydrate metabolism
    dGluc = v1 - v2 - Gluc*Vratio
    dG6P = v2 - v3 - v10 - G6P*Vratio
    dFBP = v3 - v4 - FBP*Vratio
    dGAP = v4 + v5 - v6 - GAP*Vratio
    dDHAP = v4 - v5 - v11 - DHAP*Vratio
    dPyr = v6 - v7 - v8 - v9 - Pyr*Vratio
    dG3P = v11 - v12 - G3P*Vratio
    dGlycin = v12 - v13 - Glycin*Vratio ##KS
    dGlycex = v13/kx
    dNADH = v6 + 4*v7 + v8 + v11 + v14 - v15 - NADH*Vratio
    dNADp = -v6 - 4*v7 - v8 - v11 - v14 + v15 - NADp*Vratio
    dATP = -v2 - v3 + 2*v6 + v7 - v9 - v10 + 3*v14 - v16 - ATP*Vratio
    dADP = v2 + v3 - 2*v6 - v7 + v9 + v10 - 3*v14 + v16 - ADP*Vratio
    
        # volume and osmotic pressure  
    dc = v12 - v13 - c*Vratio
    dPIi = PIi*dc/c
    dVos = -G*Lp*(PIt+PIe-PIi)
    dVratio = dVos/Vos


    # **********
    # Other dependancies
    # **********

    k13 = k13_0 * math.pow((PIt / PIt_0), n13)

    # phosphorelay module
    Sln1 = Sln1total - Sln1P  ###KS
    Ypd1 = Ypd1total - Ypd1P  ###KS
    Ssk1 = Ssk1total - Ssk1P  ###KS

    # MAP kinase cascade
    Ssk2P = Ssk2total - Ssk2
    Pbs2P2 = Pbs2total - (Pbs2 + Pbs2P)
    Hog1T = (Hog1 + Hog1P + Hog1P2cyt) * (Vcyt / Vtotal) + (Hog1P2nuc + Hog1nuc) * (Vnuc / Vtotal)  ### KS

    # transcritpion and translation
    Glyc = Glycin + Glycex

    # carbohydrate metabolism

    # volume and osmotic pressure
    dVtotal = dVos
    dVcyt = dVos
    Vos = Vtotal - Vb
    Vratio = dVos / Vos
    # Pie = Pife*(1-np.exp(-t/lamdam) + Pie_0*np.exp(-t/lambdam))  what is Pife??????
    PIt = PIt_0 * (1 - (Vos_0 - Vos) / (Vos_0 - VPIt_0))


        
    return (dc, PIt, dPIi, PIe, dVtotal, Vb, Vnuc, dVcyt, dVratio, dVos, VPIt, dSln1P, dYpd1P, dSsk1P, Sln1total, Ypd1total, Ssk1total, dSln1, dYpd1, dSsk1, dSsk2, dPbs2, dPbs2P, dHog1, dHog1P, dHog1P2cyt, dHog1P2nuc, dHog1nuc, Ssk2total, Pbs2total, Hog1total, dSsk2P, dPbs2P2, Hog1T, dmRNAnuc1, dmRNAcyt1, dProtein1, dGluc, dG6P, dFBP, dGAP, dDHAP, dPyr, dG3P, dGlycin, dGlycex, dNADH, dNADp, dATP, dADP)

# ------------------------------------        
# Solve using odeint
# ------------------------------------
solution = odeint(myODE,ini,t,args=(params,)) 

    # phosphorelay module
tSln1 = solution[:,0]
tSln1P = solution[:,1]
tYpd1 = solution[:,2]
tYpd1P = solution[:,3]
tSsk1 = solution[:,4]
tSsk1P =solution[:,5]

    # MAP kinase cascade
tSsk2 = solution[:,6]
tSsk2P = solution[:,7]
tPbs2 = solution[:,8]
tPbs2P = solution[:,9]
tPbs2P2 = solution[:,10]
tHog1 = solution[:,11]
tHog1P = solution[:,12]
tHog1P2cyt = solution[:,13]
tHog1P2nuc = solution[:,14]
tHog1nuc = solution[:,15]
    
    # transcritpion and translation
tmRNAnuc1 = solution[:,16]
tmRNAcyt1 = solution[:,17]
tProtein1 = solution[:,18]
#tmRNAnuc1 = solution[:,16]
#tmRNAcyt1 = solution[:,17]
#tProtein1 = solution[:,18]
    
    # carbohydrate metabolism
tGluc = solution[:,19]
tG6P = solution[:,20]
tFBP = solution[:,21]
tGAP = solution[:,22]
tDHAP = solution[:,23]
tPyr = solution[:,24]
tG3P = solution[:,25]
tGlycin = solution[:,26]
tGlycex = solution[:,27]
tNADH = solution[:,28]
tNAD = solution[:,29]
tATP = solution[:,30]
tADP = solution[:,31]
tProtein1_0 = solution[:,32] ##KS

# volume and osmotic pressure  

    
# ------------------------------------
# Show over time
# ------------------------------------
fig, ax = plt.subplots()

    # phosphorelay module
ax.plot(t, tSln1, label='Sln1')
ax.plot(t, tSln1P, label='Sln1P')
ax.plot(t, tYpd1, label='Ypd1')
ax.plot(t, tYpd1P, label='Ypd1P')
ax.plot(t, tSsk1, label='Ssk1')
ax.plot(t, tSsk1P, label='Ssk1P')

    # MAP kinase cascade
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
    
    # transcritpion and translation
ax.plot(t, tmRNAnuc1, 'r', label='mRNAnuc1')
ax.plot(t, tmRNAcyt1, 'g', label='mRNAcyt1')
ax.plot(t, tProtein1, 'y', label='Protein1')
    
    # carbohydrate metabolism
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

# volume and osmotic pressure  

legend = ax.legend(loc='center', shadow=True)
plt.show()