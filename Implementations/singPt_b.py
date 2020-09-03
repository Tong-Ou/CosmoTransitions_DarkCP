from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:33:49 2018

@author: yik
"""
import sys
sys.path.append('/home/tong/Chicago/EWPhT/cosmotransition_z2sb/cosmoTransitions/')

import baseMo_s_b_cwd as bm

import matplotlib.pyplot as plt

import numpy as np

from cosmoTransitions import helper_functions
  
m = bm.model(1.401946636699999999e+03, -6.978482769999999391e+01, 5.757000000000000312e-02, 3.131300000000000056e-03, -3.626290000000000069e-02 ,v2re = 1000.**2.)

print("\n")
print("\n")

print("The T=0 potential of the model reads")

print("\n")

print("Now let's find the phase transitions:")

m.calcTcTrans()

print("\n \n All the phase transitions of such a model are")

m.prettyPrintTcTrans()

print("And the T-dependent 'phase norm' reads")

plt.figure()
m.plotPhasesPhi()
plt.show()

plt.figure()
m.plotPhases2D()
plt.show()

"""
Note: to be completed: 
models may have probolems calculating tunneling (possibly due to the strength of phase transitions).
We need Tc info instead of Tn info. 
So such a step shall be neglected at this point.

"""
print("\n \n")

print("Now let's find the corresponding tunneliings:")

m.findAllTransitions()

print("\n \n All the tunnelings/phase transitions of such a model are")

m.prettyPrintTnTrans()


def physpara(m):
    
    vp = m.findMinimum(T=0.)
    
    v2phy = vp[0]**2.
    
    tanbphy = vp[1]/vp[0]
    
    m2phys = m.d2V(vp, T=0.)
        
    m2physd, eigv = np.linalg.eig(m2phys)   
    
    if all((m2physd[0] >= 0., m2physd[1] >= 0.)):
                      
        sintphy = eigv[0][1]
           
        m1phy, m2phy = m2physd**.5
    
    else:
        
        sintphy = 0.
                
        m1phy = 0.
        
        m2phy = 0.                    
        
    return [v2phy, tanbphy, m1phy, m2phy, sintphy]


def massbase(m):
    
    vp = m.findMinimum(T=0.)
    
    m2phys = m.d2V(vp, T=0.)
    
    m2physd, eigv = np.linalg.eig(m2phys)
    
    vpm = np.matmul(eigv.T, vp)
    
    return eigv, vpm


def Vtotm(X, T, m):
    
    T = np.asanyarray(T, dtype=float)
    X = np.asanyarray(X, dtype=float)

    hm,sm = X[...,0], X[...,1]
    
    eigv, vpm = massbase(m)
    
    hg = eigv[0][0]*(hm + vpm[0]) + eigv[0][1]*(sm + vpm[1])
    
    sg = eigv[1][0]*(hm + vpm[0]) + eigv[1][1]*(sm + vpm[1])
    
    xg = np.array([hg, sg])
    
    xg = np.rollaxis(xg, 0, len(xg.shape))
            
    return m.Vtot(xg, T)

        
def Dh(X,T, m):
    return helper_functions.gradientFunction(
                Vtotm, m.x_eps, m.Ndim, m.deriv_order)(X,T, m)[...,0]


def Ds(X,T, m):
    return helper_functions.gradientFunction(
                Vtotm, m.x_eps, m.Ndim, m.deriv_order)(X,T, m)[...,1]


def Dhh(X,T, m):
    return helper_functions.gradientFunction(
                Dh, m.x_eps, m.Ndim, m.deriv_order)(X,T, m)[...,0]


def Dhs(X,T, m):
    return helper_functions.gradientFunction(
                Dh, m.x_eps, m.Ndim, m.deriv_order)(X,T, m)[...,1]


def Dhss(X,T, m):
    return helper_functions.gradientFunction(
                Dhs, m.x_eps, m.Ndim, m.deriv_order)(X,T, m)[...,1]
    
def nst(x):
    
    nst = 0.
    
    muh2, mus2, lh, ls, lm = x[0], x[1], x[2], x[3], x[4]
    
    m = bm.model(muh2, mus2, lh, ls, lm, v2re = 1000.**2.)
    
    vphys = physpara(m)[0]**.5
    
    tanbphy = physpara(m)[1]
            
    m1phys = physpara(m)[2]
            
    m2phys = physpara(m)[3]
    
    sintphy = physpara(m)[4]
            
    print '.',
                                        
    if all((vphys <= 248., vphys >= 244., tanbphy >= 0.001)):
                
        if all((m1phys <= 127., m1phys >= 123., sintphy <= .4)) or all((m2phys <= 127., m2phys >= 123., (1. - sintphy**2.)**.5 <= .4)):
            
            m.calcTcTrans()
    
            trans = m.TcTrans
    
    
            for k in range(len(trans)):
                        
                tc = trans[k]['Tcrit']
                sh = abs(trans[k]['low_vev'][0]-trans[k]['high_vev'][0])
                        
                if trans[k]['trantype'] == 1 and -sh/tc <= nst and trans[k]['high_vev'][0] <= 2.:
                        
                    nst = -sh/tc
    
    if nst <= -0.9:
        
        fnst = -0.9
    
    else:
        
        fnst = nst
    
                    
    print nst
            
    return fnst, nst



for i in nd:
    m = bm.model(i[3],i[4],i[0],i[1],i[2],v2re = 1000.**2.)
    imax = optimize.fmin(f,[500., 1000.],disp=0)
    ptmax = np.append(i, imax, axis=0)
    ndmax.append(ptmax)
                
    
