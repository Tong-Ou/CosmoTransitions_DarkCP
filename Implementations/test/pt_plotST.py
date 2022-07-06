#from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:33:49 2018

@author: Tong Ou (adapted from Yikun's codes)
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

t0 = time.clock()

ntasks = int(sys.argv[1])
FILE = sys.argv[2]
OUT_PATH = sys.argv[3]

if not os.path.isdir(OUT_PATH):
    os.makedirs(OUT_PATH)

filename = '%s.npy' % FILE
paras = np.load(filename, allow_pickle = True)

#fig = plt.figure()
'''
scan_task = range(len(paras))
#scan_task = [40, 1, 63, 44, 25, 47, 57, 17]
rank_task = scan_task[rank:len(scan_task):numtasks]
logfile = '%s/pt_%s.log' % (OUT_PATH, rank)
log = open(logfile, 'w')
sys.stdout = log
phase_dict = {}
'''
#rank_task = [210, 19, 468, 148, 336, 276, 127, 287]
rank_task = [8]
for index in rank_task:

    para = paras[index]
    print('The parameters are:')
    print( 'Index: %s vh2:%s vs2:%s lh:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s v2re:%s' % (index, para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9]))
    
    mt = bmt.model(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9])
    #print (mt.Vtot([100.,100.,100.], T=0.))
    
    vphy = para[10]
    zeroTLocalMin = para[12]

    findex = index % ntasks
    phases = dd.io.load('%s/phases_%s.h5' % (OUT_PATH, findex))[index]

    '''    
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
    phase_dict.update({index:phases})
    
    mt.calcTcTrans(vphy, zeroTLocalMin)
    
    print("\n \n All the phase transitions of such a model are")
    
    mt.prettyPrintTcTrans()

    print('The T-dependent phase potential reads')
    mt.plotPhasesV()
    plt.savefig('%s/V_T_%s.pdf' % (OUT_PATH, index))
    plt.clf()
    
    print("And the T-dependent 'phase norm' reads")
    
    #plt.figure()
    mt.plotPhasesPhi()
    #plt.show()
    plt.savefig('%s/phi_T_%s.pdf' % (OUT_PATH, index))
    plt.clf()
    
    #plt.figure()
    mt.plotPhases2D()
    #plt.show()
    plt.savefig('%s/vs_vh_%s.pdf' % (OUT_PATH, index))
    plt.clf()

    mt.plotPhase2DS()
    plt.savefig('%s/vs_va_%s.pdf' % (OUT_PATH, index))
    plt.clf()
    
    """
    Note: to be completed: 
    models may have probolems calculating tunneling (possibly due to the strength of phase transitions).
    We need Tc info instead of Tn info. 
    So such a step shall be neglected at this point.
    """

    print("\n \n")
    
    print("Now let's find the corresponding tunneliings:")
    
    try:
        mt.findAllTransitions(vphy, zeroTLocalMin)
    except KeyError as err:
        print ('Skipping due to KeyError: %s...' % err)
        continue
    except ValueError as err:
        print ('Skipping due to ValueError: %s...' % err)
        continue
    except:
        print ('Skipping due to unexpected error...')
        continue
    
    print("\n \n All the tunnelings/phase transitions of such a model are")

    ses, V, dV, makePlot=False,                                                                                                       phitol=1e-8, overlapAngle=45.0,                                                                                                        nuclCriterion=lambda S,T: S/(T+1e-100) - 140.0,                                                                                        verbose=True, OUT_PATH='.', index=0,                                                                                                   fullTunneling_params={})mt.prettyPrintTnTrans()
    '''

    mt.plotNuclCriterion(phases, OUT_PATH, index)
    #plt.savefig('%s/S_T_%s.png' % (OUT_PATH, index))
    #plt.clf()
    

# Save phases
#dd.io.save('%s/phases_%s.h5' % (OUT_PATH, rank), phase_dict)

t1 = time.clock()
print ('Run time: %s' % (t1-t0))
