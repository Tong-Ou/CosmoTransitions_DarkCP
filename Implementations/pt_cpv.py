from __future__ import division

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 20:33:49 2018

@author: yik
"""

import sys
sys.path.append('/home/tong/Work/EWPhT/cosmotransition_z2s/cosmoTransitions/')

import os
#import baseMo_s_t as bmt
import baseMo_s_cpv as bmt

import matplotlib.pyplot as plt
import numpy as np

import mpi4py.MPI as mpi

comm = mpi.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

FILE = sys.argv[1]
OUT_PATH = sys.argv[2]

#paras_0 = np.load('%s_0.npy' % FILE)
paras_1 = np.load('%s_1.npy' % FILE, allow_pickle = True) #unsafe
paras_2 = np.load('%s_2.npy' % FILE, allow_pickle = True)
#paras_3 = np.load('%s_3.npy' % FILE)
paras = np.concatenate((paras_1, paras_2),axis = 0)
np.save('%s.npy' % FILE, paras)

fig = plt.figure()

for index in range(len(paras)):
    if index%size != rank:
        continue
    else:
        logfile = '%s/pt_%s.log' % (OUT_PATH, rank)
        if index==rank and os.path.exists(logfile):
           os.remove(logfile)
        log = open(logfile, 'a')
        sys.stdout = log

        para = paras[index]
        print('The parameters are:')
        print( 'Index: %s vh2:%s vs2:%s lh:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s v2re:%s' % (index, para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9]))
        
        #mt = bmt.model(0.3,0.036,-0.171, 125., 246.**2., 1000.**2.)
        mt = bmt.model(para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9])
        vphy = np.array([para[10], 0., 0.])
        zeroTLocalMin = para[12]
        
        print("\n")
        print("\n")
        
        print("The T=0 potential of the model reads")
        
        bmt.vsh(mt, [-300., 300., -400., 400., 0.], 0., cmap='RdGy')
        plt.savefig('%s/V0_%s.png' % (OUT_PATH, index))
        plt.clf()
        
        print("\n")
        print("\n")
         
        print("Now let's find the phase transitions:")
        
        mt.calcTcTrans(vphy, zeroTLocalMin)
        
        print("\n \n All the phase transitions of such a model are")
        
        mt.prettyPrintTcTrans()
        
        print("And the T-dependent 'phase norm' reads")
        
        #plt.figure()
        mt.plotPhasesPhi()
        #plt.show()
        plt.savefig('%s/phi_T_%s.png' % (OUT_PATH, index))
        plt.clf()
        
        #plt.figure()
        mt.plotPhases2D()
        #plt.show()
        plt.savefig('%s/vs_vh_%s.png' % (OUT_PATH, index))
        plt.clf()

        mt.plotPhase2DS()
        plt.savefig('%s/vs_va_%s.png' % (OUT_PATH, index))
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
        
        '''
        mt.plotNuclCriterion()
        plt.savefig('%s/S_T_%s.png' % (OUT_PATH, index))
        plt.clf()
        '''
