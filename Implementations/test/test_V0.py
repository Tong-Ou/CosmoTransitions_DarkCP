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

#import baseMo_s_t as bmt
import baseMo_s_cpv as bmt

import numpy as np

#import mpi4py.MPI as mpi

#comm = mpi.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
'''
FILE = sys.argv[1]

if not os.path.isdir(OUT_PATH):
    os.mkdir(OUT_PATH)

filename = '%s.npy' % FILE
paras = None
if os.path.exists(filename):
    paras = np.load(filename, allow_pickle = True)
else:
    ncpu = 100
    para_list = []
    for i in range(ncpu):
        filename = '%s_%s.npy' % (FILE, i)
        if os.path.exists(filename):
            para = np.load(filename, allow_pickle = True)
            para_list.append(para)
        else:
            continue
    paras = np.concatenate(para_list,axis = 0)
    np.save('%s.npy' % FILE, paras)
'''
#fig = plt.figure()
'''
scan_task = range(len(paras))
rank_task = scan_task[rank:len(paras):numtasks]
logfile = '%s/pt_%s.log' % (OUT_PATH, rank)
log = open(logfile, 'w')
sys.stdout = log
'''
for index in range(1):

    para = [60516, 100.**2, 0.129, 1.3343720573, 9.8347897898, 3000.0, 2.0, 4.6058404343, 100., 600.**2]
    print('The parameters are:')
    print( 'Index: %s vh2:%s vs2:%s lh:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s v2re:%s' % (index, para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9]))
    
    mt = bmt.model(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9])
    print ('Vtot(vh):')
    print (mt.Vtot([246.0000260616804, -1.1415885184493248e-05, -2.246559391093995e-05], T=0.))
    print ('V0(vh):')
    print (mt.V0([246.0000260616804, -1.1415885184493248e-05, -2.246559391093995e-05]))
    print ('VCW(vh):')
    print (mt.V1([246.0000260616804, -1.1415885184493248e-05, -2.246559391093995e-05], T=0.))
    print ('Counterterm(vh):')
    print (mt.counterterm([246.0000260616804, -1.1415885184493248e-05, -2.246559391093995e-05]))

    print ('\n')
    lmin = [0, 0, 0]
    print ('Vtot(%s):' % lmin)
    print (mt.Vtot(lmin, T=0.))
    print ('V0(%s):' % lmin)
    print (mt.V0(lmin))
    print ('VCW(%s):' % lmin)
    print (mt.V1(lmin, T=0.))
    print ('Counterterm(%s):' % lmin)
    print (mt.counterterm(lmin))

    '''
    print ('Boson masses:')
    print (mt.boson_massSq([100.,100.,100.], T=0.))
    print ('Fermion masses:')
    print (mt.fermion_massSq([100.,100.,100.]))
    '''
    '''
    vphy = para[10]
    zeroTLocalMin = para[12]
    
    print("\n")
    print("\n")
    
    print("The T=0 potential of the model reads")
    
    bmt.vsh(mt, [-300., 300., -400., 400., 0.], 0., cmap='RdGy')
    plt.savefig('%s/V0_%s.pdf' % (OUT_PATH, index))
    plt.clf()
    
    print("\n")
    print("\n")
     
    print("Now let's find the phase transitions:")
    
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

    mt.prettyPrintTnTrans()
    
    
    
    mt.plotNuclCriterion()
    plt.savefig('%s/S_T_%s.png' % (OUT_PATH, index))
    plt.clf()
   '''
