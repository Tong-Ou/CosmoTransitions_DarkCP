#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 17:10:18 2018

@author: yik
"""

"""
In this module, we create a scaning over the parameter space (ms, sint, tanb) which enters the model 
defined and used in the module baseMo.py

We create a dictionary containing the parameter info and the 1st order phase transition info, including 
critical temperature, high low T vevs, and the phase transition strength. 

Usage functions are created to print/plot the scaning. 
"""

import find_n as find

import random as ran

import numpy as np

msr = [0.1, 600.]

tanbr = [0.1,10.]

sintr = [-0.3,0.3]

npt = 10000

ranscanslved = []

ranscanno = []


def getscan(msbox, tanbbox, sintbox, npt):
    msmin,msmax = msbox
    tmin, tmax = tanbbox
    smin, smax = sintbox
    for i in range(0, npt):
        ms0 = ran.uniform(msmin,msmax)
        tanb0 = ran.uniform(tmin, tmax)
        sint0 = ran.uniform(smin, smax)
        if ms0 >= 246.:
            v2re = ms0**2.
        else:
            v2re = 246.**2.
            
        try:
            
            f = find.findmhv(ms0, tanb0, sint0, v2re)
            x = f.find()
            
            if x[0] == None:
                
                ranscanno.append([ms0, tanb0, sint0])
                
                np.savetxt("outputs/ran_scan_mhv_no_n", ranscanno)
            
            else:       
            
                vp, mhp, msp = f.mhv(x[0],x[1])
                            
                ranscanslved.append([ms0, tanb0, sint0, x[0],x[1],msp,vp[1]/vp[0]])
             
                np.savetxt("outputs/ran_scan_mhv_slved_n", ranscanslved)
                
        except Exception: 
            
            print 'something wrong'
            
            pass
            
    return ranscanslved, ranscanno
                                    
scan = getscan(msr, tanbr, sintr, npt)

# np.load("ran_scan_1.npz")

#npfile.files

#npfile['check']

# npfile['scan_1pht'].item().get('sint')

"""
 
for i in allpt:
    ms0 = i[0]
    tanb0 = i[1]
    sint0 = i[2]
    v0 = i[3]
    mh0 = i[4]
    m = bm.model(ms0, tanb0, sint0, mh0, v0**2., v2re = 246.**2.)
    vp, mhp, msp = mhv(m)    
    physm = m.d2V(vp, T=0.)
    m2hs = physm[0][1]
    try:
        sintphy = sin(.5*asin((2.*m2hs)/(msp**2. - mhp**2.)))
        allptnew.append(np.append(i, sintphy))
    except:
        print (2.*physm[0][1])/(i[5]**2.-125.**2.)
        print i
"""
