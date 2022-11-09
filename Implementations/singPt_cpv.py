#from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:33:49 2018

@author: Tong Ou (modified from Yikun's codes)
"""

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

import sys

import os
import time

#import baseMo_s_t as bmt
import baseMo_s_cpv as bmt

import numpy as np
import deepdish as dd

#import mpi4py.MPI as mpi

#comm = mpi.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()

FILE = sys.argv[1]
OUT_PATH = sys.argv[2]

logfile = '%s/pt.log' % (OUT_PATH)
log = open(logfile, 'w')
sys.stdout = log

#if not os.path.isdir(OUT_PATH):
 #   os.makedirs(OUT_PATH)

filename = '%s.npy' % FILE
paras = np.load(filename, allow_pickle = True)
vphys = np.load(FILE + '_vphy.npy', allow_pickle = True)
index = 960
para = paras[index]
print('The parameters are:')
print( 'Index: %s vh2:%s vs2:%s lh:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s v2re:%s' % (index, para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9]))
    
mt = bmt.model(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9])
#print (mt.Vtot([100.,100.,100.], T=0.))
    
vphy = vphys[index]
zeroTLocalMin = np.array([[246., 0., 0.], [0., 0., abs(mt.vs2+2.*mt.ks2/mt.ls)**0.5], [1e-5, 1e-5, 1e-5]])
    
print("\n")
print("\n")
    
print("The T=0 potential of the model reads")
    
bmt.vsh(mt, [-300., 300., -400., 400., 0.], 0., cmap='RdGy')
plt.savefig('%s/V0_%s.pdf' % (OUT_PATH, index))
plt.clf()
    
print("\n")
print("\n")
     
print("Now let's find the phase transitions:")

phases = mt.getPhases(vphy, zeroTLocalMin)
    
mt.calcTcTrans(vphy, zeroTLocalMin)
    
print("\n \n All the phase transitions of such a model are")
    
mt.prettyPrintTcTrans()

print('The T-dependent phase potential reads')
mt.plotPhasesV()
plt.savefig('%s/V_T_%s.pdf' % (OUT_PATH, index))
plt.clf()
    
print("And the T-dependent 'phase norm' reads")
    
mt.plotPhasesPhi()
plt.savefig('%s/phi_T_%s.pdf' % (OUT_PATH, index))
plt.clf()
    
mt.plotPhases2D()
plt.savefig('%s/vs_vh_%s.pdf' % (OUT_PATH, index))
plt.clf()

mt.plotPhase2DS()
plt.savefig('%s/vs_va_%s.pdf' % (OUT_PATH, index))
plt.clf()
    
print("\n \n")
    
print("Now let's find the corresponding tunneliings:")
''' 
try:
    mt.findAllTransitions(vphy, zeroTLocalMin)
except KeyError as err:
    print ('Skipping due to KeyError: %s...' % err)
    pass
except ValueError as err:
    print ('Skipping due to ValueError: %s...' % err)
    pass
except:
    print ('Skipping due to unexpected error...')
    pass
    
mt.dSdT() # Compute dS/dT at Tnuc (for calculation of GW)
print("\n \n All the tunnelings/phase transitions of such a model are")
mt.prettyPrintTnTrans()
'''    
mt.plotNuclCriterion(phases, OUT_PATH, index)

