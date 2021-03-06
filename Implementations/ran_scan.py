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

import sys
sys.path.append('/home/tong/Work/EWPhT/cosmotransition_z2s/cosmoTransitions/')

#import baseMo_s_b_cw as bmcw

#import baseMo_s_b_cwd as bmcwd
import baseMo_s as bm

import deepdish as dd

import random as ran

import numpy as np

import time

import os


import mpi4py.MPI as mpi

import math

import matplotlib.pyplot as plt

import time

comm = mpi.COMM_WORLD

FILE = sys.argv[1]

t0 = time.clock()

#lhr = [0.01,0.05]
#For tree-level potential, lh=mh^2/2vH^2
lh = 0.129
vh = 246.

lsr = [0., 4*np.pi]

lshr = [0., 4*np.pi]

ydr = [1e-2, 10.]

#vh2logr = [0.,4.]

#vs2logr = [0.,4.]
#vh2 = [230,260]

vsr = [100., 500.]

#cHr = [0.,1.]

#csr = [0.,1.]

#m12r = [-100000., 100000.]

#m22r = [-10000., 0.]

npt = int(sys.argv[2])

scan = []

BASE_PATH = '/home/tong/Work/EWPhT/cosmotransition_z2s/Implementations/'

def physpara(m):
    
    vp = m.findMinimum(T=0.)
    print ('Minimum: [%s, %s]' % (vp[0],vp[1]))

    #check whether this is global minimum
    #VS = np.sqrt(m.vs2)
    #vs = m.findMinimum([np.array([0.1, VS])], T=0.)
    #if m.Vtot(vs, T=0.) < m.Vtot(vp, T=0.):
     #   return None
    
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
        
        

def trans(i, m):

    print "\n check point %s \n" %[i, m.lh, m.ls, m.lsh, m.vh2, m.vs2] 
    
    m.calcTcTrans() 

    trans = m.TcTrans

    sfo = False

    fo = False
    
    check = False

    foph = []
    
    sfoph = []
    
    for k in range(len(trans)):
        tc = trans[k]['Tcrit']
        sh = abs(trans[k]['low_vev'][0]-trans[k]['high_vev'][0])/tc
        if trans[k]['trantype'] ==1:
            foph.append(sh)
            fo = True
            if sh >= 1.: 
                sfoph.append([tc, sh])
                sfo = True
               
    for key in m.phases:
        if m.phases[key].check:               
            check = True
    
    #return trans, sfoph, check
    return fo, foph, sfo, sfoph
    
def t0vev(m):
    
    wvev = m.Vtot(m.findMinimum(T=0.),T=0.) - m.Vtot(m.findMinimum(X=[0.,0.],T=0.),T=0.) > 1.
    
    return wvev


def thvev(m):
    
    htX =  m.findMinimum(T=1000.)
    
    wvev = (abs(htX[...,0]) > 10.**10.) or (abs(htX[...,1]) > 10.**10.)
    
    return wvev


def getscani(i, m):
    
    
    if any([m.lh > 4.*np.pi/3., m.ls > 4.*np.pi/3., m.lsh > 16.*np.pi/3.]): 
           
        print('wrong paras')
        
        scani = None
        
        
    else: 
                    
        scani = []
                         
        scani.append([m.lh, m.ls, m.lsh, m.vh2, m.vs2])
        
        scani.append(physpara(m))
            
        if t0vev(m) or thvev(m):
            
            scani.append([])
            
            scani.append([])
            
            scani.append(True)
                
        else:

            transit, sfoph, check = trans(i, m)
            
            scani.append([sfoph])
            
            scani.append([transit])
            
            scani.append(check)
         
    return scani
                
                                                  

def getscan(l2box, lmbox, m2box, ydbox, npt):
 
   # l1min,l1max = l1box
    l2min, l2max = l2box
    lmmin, lmmax = lmbox
    #m12min, m12max = m12box
    m2min, m2max = m2box
    ydmin, ydmax = ydbox
    #chmin, chmax = chbox
    #csmin, csmax = csbox
        
    scan_task = range(npt)
    
    rank = comm.Get_rank()

    size = comm.Get_size()
        
    scan_rank = []
                
    ran.seed(time.time() + rank)
    #plt.figure(figsize=(12,5))

    ls_scan = []
    lsh_scan = []
    yd_scan = []
    ls_sfo = []
    lsh_sfo = []
    yd_sfo = []
    ls_fo = []
    lsh_fo = []
    yd_fo = []
    #sh_fo = []
        
    for n in scan_task:
        
        if n%size != rank: 
           continue
       
        else:
            logfile = '%s/%s_%s.log' % (BASE_PATH, FILE, rank)
            if n==rank and os.path.exists(logfile):
                os.remove(logfile)
            log = open(logfile, 'a')
            sys.stdout = log

            # Tong: generate random number between the given min and max
            # Why not create a np.linspace and scan one by one? - Not realistic for 5 parameters
            #l1 = ran.uniform(l1min,l1max)
            l1 = lh
            l2 = ran.uniform(l2min,l2max)
            #l2 = ran.uniform(0, 5.)
            lm = ran.uniform(lmmin,lmmax)
            #m1log = ran.uniform(m1min,m1max)
            m2log = ran.uniform(m2min,m2max)
            #m2log = ran.uniform(200, 250)
            yd = ran.uniform(ydmin, ydmax)

            m12 = vh**2
            m22 = m2log**2
            #yd = np.sqrt(yd2)
            #ch = ran.uniform(chmin,chmax)
            #cs = ran.uniform(csmin,csmax)
            #
            #m12 = 10.**m12log
            #
            #m22 = 10.**m22log
                      
            v2re = 600.**2.
    
            # Zero-T global minimum condition
            if l2*m22**2 >= l1*m12**2:
                #scan_task += 1
                continue
                
            mcwd = bm.model(m12, m22, l1, l2, lm, yd, v2re)

            # Physical parameters
            phy = physpara(mcwd)
            if phy is None:
                continue
            vphy, tanbphy, m1phy, m2phy, sintphy = phy[0]**.5, phy[1], phy[2], phy[3], phy[4]

            #Tong: Constraints: 1. Higgs mass and vev.
                        
            if all((vphy <= 248., vphy >= 244.)):
                
                if all ((m1phy <= 127., m1phy >= 123.)):

                    # 2. FO condition
                    fo, foph, sfo, sfoph = trans(n, mcwd)
               
                #if all((m1phy <= 127., m1phy >= 123., sintphy <= .4)) or all((m2phy <= 127., m2phy >= 123., (1. - sintphy**2.)**.5 <= .4)):
                                    
                    if sfo:
                        ch = mcwd.ch()
                        cs = mcwd.cs()
                        if cs*l1*m12 >= ch*l2*m22:
                            diff1 = cs*l1*m12 - ch*l2*m22
                            print ('\n')
                            print ('cs*lh*vh2 - ch*ls*vs2 = %s' % diff1)
                            print ('WARNING: Strong first-order transition found without satisfying high-T conditions.')
                        else:
                            print ('Find strong first-order transition!!')
                        """
                        scan_rank.append([m12, m22, l1, l2, lm, yd, v2re, vphy, tanbphy, m1phy, m2phy, sintphy])
              
                        filename = '%s_%s' % (FILE, rank)
              
                        np.save(os.path.join(BASE_PATH, filename), scan_rank)
                        """
                        ls_sfo.append(l2)
                        lsh_sfo.append(lm)
                        yd_sfo.append(yd)

                    if fo:
                        if not sfo:
                            yd_fo.append(yd)
                            # v/Tc
                            # sh_fo.append(foph[0])
                            ls_fo.append(l2)
                            lsh_fo.append(lm)
                    
                        scan_rank.append([m12, m22, l1, l2, lm, yd, v2re, vphy, tanbphy, m1phy, m2phy, sintphy])
                        filename = '%s_%s' % (FILE, rank)
                        np.save(os.path.join(BASE_PATH, filename), scan_rank)

                    else:
            
                        ls_scan.append(l2)
                        lsh_scan.append(lm)
                        yd_scan.append(yd)
                    
                else:
                                        
                    pass
           
            else:
                
                pass

    para_dict = {'ls_scan':ls_scan, 'lsh_scan':lsh_scan, 'yd_scan':yd_scan,
                 'ls_fo':ls_fo, 'lsh_fo':lsh_fo, 'yd_fo':yd_fo,
                 'ls_sfo':ls_sfo, 'lsh_sfo':lsh_sfo, 'yd_sfo':yd_sfo}
    filename = '%s_%s_plot.h5' % (FILE, rank)
    dd.io.save(os.path.join(BASE_PATH, filename), para_dict)

                                                 
scan = getscan(lsr, lshr, vsr, ydr, npt)
tf = time.clock()
dt = tf-t0
print ('Run time: %s' % dt)


