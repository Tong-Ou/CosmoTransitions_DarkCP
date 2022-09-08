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
with open('./nuc_cpvparas.dat', 'r') as f:
    lines = f.readlines()

rank = int(sys.argv[1]) #task_id
v2re = 1e6

numtasks = int(sys.argv[2])

#FILE = sys.argv[3]
OUT_PATH = sys.argv[3]
scan_task = range(len(lines))
rank_task = scan_task[rank:len(scan_task):numtasks]
logfile = '%s/nuc_%s.log' % (OUT_PATH, rank)
log = open(logfile, 'w')
sys.stdout = log

for i in rank_task:

    paras = lines[i].split()
    index, vs2, ls, lsh, ks2, m0, yd, theta = int(paras[0]), float(paras[1]), float(paras[2]), float(paras[3]), float(paras[4]), float(paras[5]), float(paras[6]), float(paras[7])
    print('The parameters are:')
    print( 'Index: %s vh2:%s vs2:%s lh:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s v2re:%s' % (index, 246**2,vs2,0.129,ls,lsh,ks2,yd,theta,m0,v2re))

    if ls > 4*np.pi:
        print ('ls > 4 Pi, skipping...')
        continue
    mt = bmt.model(246.**2, vs2, 0.129, ls, lsh, ks2, yd, theta, m0, v2re)
    #print (mt.Vtot([100.,100.,100.], T=0.))

    vphy = np.array([246., 0., 0.])
    zeroTLocalMin = np.array([[246., 0., 0.], [0., 0., abs(mt.vs2+2.*mt.ks2/mt.ls)**0.5], [1e-5, 1e-5, 1e-5]])
    #zeroTLocalMin = lmins[index]

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
    #phase_dict.update({index:phases})

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

    mt.dSdT() # Compute dS/dT at Tnuc (for calculation of GW)
    print("\n \n All the tunnelings/phase transitions of such a model are")
    mt.prettyPrintTnTrans()

t1 = time.clock()
print ('Run time: %s' % (t1-t0))
