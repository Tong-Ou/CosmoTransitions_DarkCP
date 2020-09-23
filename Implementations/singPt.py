from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:33:49 2018

@author: yik
"""

import sys
sys.path.append('/home/tong/Chicago/EWPhT/cosmotransition_z2s/cosmoTransitions/')

#import baseMo_s_t as bmt
import baseMo_s as bmt

import matplotlib.pyplot as plt
import numpy as np

paras = np.load('outputs/test.npy')
para = paras[33]
print('The parameters are:')
print 'vh2:%s vs2:%s lh:%s ls:%s lsh:%s yd:%s v2re:%s' % (para[0],para[1],para[2],para[3],para[4],para[5],para[6])

#mt = bmt.model(0.3,0.036,-0.171, 125., 246.**2., 1000.**2.)
mt = bmt.model(para[0],para[1],para[2],para[3],para[4],para[5],para[6])

print("\n")
print("\n")

print("The T=0 potential of the model reads")

bmt.vsh(mt, [-300., 300., -400., 400.], 0., cmap='RdGy')

print("\n")
print("\n")

print("Now let's find the phase transitions:")

mt.calcTcTrans()

print("\n \n All the phase transitions of such a model are")

mt.prettyPrintTcTrans()

print("And the T-dependent 'phase norm' reads")

plt.figure()
mt.plotPhasesPhi()
#plt.show()

plt.figure()
mt.plotPhases2D()
#plt.show()

"""
Note: to be completed: 
models may have probolems calculating tunneling (possibly due to the strength of phase transitions).
We need Tc info instead of Tn info. 
So such a step shall be neglected at this point.
"""
print("\n \n")

print("Now let's find the corresponding tunneliings:")

mt.findAllTransitions(makePlot=True)

print("\n \n All the tunnelings/phase transitions of such a model are")

mt.prettyPrintTnTrans()

plt.figure()
mt.plotNuclCriterion()
plt.show()


