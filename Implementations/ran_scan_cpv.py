#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 17:10:18 2018

@author: Tong Ou (adapted from Yikun's codes)
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
import baseMo_s_cpv as bm

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

lsr = [0., 6.]

lshr = [0., 4*np.pi]

ks2r = [0., 100.]

#ydr = [1e-2, 10.]
ydr = [1e-2, 10.]
thetaYr = [0., 2*np.pi]

#m0r = [1e-3, 1e3]
m0r = [1e-3, 100.]

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
    print ('\n')
    print ('ls:%s lsh:%s vs2:%s ks2:%s m0:%s yd:%s thetaYd:%s ' % (m.ls, m.lsh, m.vs2, m.ks2, m.m0, m.yd, m.thetaY))
    print ('non-tachyonic: %s' % bool((0.5*m.lsh*m.vh2 - 2*m.ks2)/m.vs2-m.ls>0.))
    print ('condition 2-a: %s' % all(((m.lh*m.vh2**2-4.*m.ks2*m.vs2+np.sqrt(m.lh**2*m.vh2**4-8.*m.lh*m.vh2**2*m.ks2*m.vs2))/(2.*m.vs2**2)-m.ls>0., (m.lh*m.vh2**2-4.*m.ks2*m.vs2-np.sqrt(m.lh**2*m.vh2**4-8.*m.lh*m.vh2**2*m.ks2*m.vs2))/(2.*m.vs2**2)-m.ls<0. )))

    '''
    gradV_exp = m.gradV([246.,0.,0.], T=0.)
    d2V_exp = m.d2V([246.,0.,0.], T=0.)
    m2sq_exp, eig_exp = np.linalg.eig(d2V_exp)
    print ('\n')
    print ('First derivative at expected Higgs minimum: %s' % gradV_exp)
    print ('Mass eigenvalues at expected Higgs minimum: %s' % m2sq_exp)
    if any((m2sq_exp[0]<0, m2sq_exp[1]<0, m2sq_exp[2]<0)):
        print ('Expected Higgs vev is a saddle point, exiting...')
        return None
    '''
    vp = m.findMinimum(T=0.)
    print ('Higgs minimum: [%s, %s, %s]' % (vp[0],vp[1],vp[2]))
    if vp[0] > 248. or vp[0] < 244.:
        print ('Not at the expected Higgs vev, exiting...')
        return 1
    else:
        m2phys = m.d2V(vp, T=0.)
        #mh2 = -m.lh*m.vh2+3.*m.lh*vp[0]**2.+0.5*m.lsh*(vp[1]**2.+vp[2]**2.)
        #ms2 = -m.ls*m.vs2+3.*m.ls*vp[1]**2.+m.ls*vp[2]**2.+0.5*m.lsh*vp[0]**2.+2.*m.ks2
        #ma2 = -m.ls*m.vs2+3.*m.ls*vp[2]**2.+m.ls*vp[1]**2.+0.5*m.lsh*vp[0]**2.-2.*m.ks2
        #msh2 = m.lsh*vp[0]*vp[1]
        #msa2 = 2.*m.ls*vp[1]*vp[2]
        #mah2 = m.lsh*vp[0]*vp[2]
        #m2phys = np.array([[mh2, msh2, mah2],[msh2, ms2, msa2],[mah2, msa2, ma2]], dtype=float)
        m2physd, eigv = np.linalg.eig(m2phys)
        # Double check
        if any((m2physd[0] < 0., m2physd[1] < 0., m2physd[2] < 0.)):
            print ('Returned Higgs vacuum is a saddle point, exiting...')
            return 1
        else:
            pass
        m1phy = m2physd[0]**0.5
        if m1phy < 123. or m1phy > 127.:
            print ('Higgs mass %s is not as expected, exiting...' % m1phy)
            return 2
        else:
            print('Search for other local minima and check whether this is global minimum...')
            vh = 246.
            vs = (m.vs2+2.*m.ks2/m.ls)**0.5
            vscan = [0.1, vh, vs, 1e5]
            minX = []
            for v1 in vscan:
                for v2 in vscan:
                    for v3 in vscan:
                        lmin = m.findMinimum(np.array([v1, v2, v3]), T=0.)
                        if m.Vtot(lmin, T=0.) < m.Vtot(vp, T=0.):
                            if np.sum((abs(lmin)-abs(vp))**2 ,-1)**.5 >= 0.01:
                                print ('Local minimum %s is lower than Higgs minimum' % lmin)
                                print ('Higgs minimum is not global minimum, exiting...')
                                return 3
                            else:
                                pass
                        else:
                            pass
                        minX.append(lmin)
   
            
            '''
            tanbphy = vp[1]/vp[0]
            
            if any((m1phy == None, eigv == None)):

                m1phy = 0.

                sintphy = 0.
            
            else:

                sintphy = eigv[0][1]
            '''
                
            return [vp, m1phy, minX]
        
        
def trans(i, m):

    print "\n check point %s \n" %[i, m.lh, m.ls, m.lsh, m.ks2, m.vh2, m.vs2] 
    
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

'''
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
                
 '''                                                 

def getscan(l2box, lmbox, ks2box, vsbox, ydbox, thetaYbox, m0box, npt):
 
   # l1min,l1max = l1box
    l2min, l2max = l2box
    lmmin, lmmax = lmbox
    ks2min, ks2max = ks2box
    #m12min, m12max = m12box
    vsmin, vsmax = vsbox
    ydmin, ydmax = ydbox
    thetaYmin, thetaYmax = thetaYbox
    m0min, m0max = m0box
    #chmin, chmax = chbox
    #csmin, csmax = csbox
        
    scan_task = range(npt)
    
    rank = comm.Get_rank()

    size = comm.Get_size()
        
    scan_rank = []
                
    ran.seed(time.time() + rank)
    #plt.figure(figsize=(12,5))

    ls_sfo = []
    lsh_sfo = []
    yd_sfo = []
    ks2_sfo = []
    ls_fo = []
    lsh_fo = []
    yd_fo = []
    ks2_fo = []
    #sh_fo = []\
    ls_novh = []
    lsh_novh = []
    yd_novh = []
    m0_novh = []
    ls_nomh = []
    lsh_nomh = []
    yd_nomh = []
    m0_nomh = []
    ls_vs = []
    lsh_vs = []
    yd_vs = []
    m0_vs = []
    ls_vh = []
    lsh_vh = []
    yd_vh = []
    m0_vh = []
        
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
            #ks2 = ran.uniform(ks2min, ks2max)
            ks2 = 3000.
            #m1log = ran.uniform(m1min,m1max)
            #vslog = ran.uniform(vsmin,vsmax)
            vslog = 100.
            #m2log = ran.uniform(200, 250)
            yd = ran.uniform(ydmin, ydmax)
            thetaY = ran.uniform(thetaYmin, thetaYmax)
            m0 = ran.uniform(m0min, m0max)

            vh2 = vh**2
            vs2 = vslog**2
                     
            v2re = 600.**2.
    
            # Zero-T global minimum condition (valid at tree-level)
            '''
            if l2*m22**2 >= (l1*m12**2-2*ks2*m22-ks2**2/l2):
                #scan_task += 1
                continue                
               ''' 

            # Zero-T non-tachyonic conditions at tree-level
            '''
            print ('\n')
            if any((-l2*vs2 + 0.5*lm*vh2 + 2.*ks2 < 0, -l2*vs2 + 0.5*lm*vh2 - 2.*ks2 < 0)):
                ls_scan.append(l2)
                lsh_scan.append(lm)
                print ('Not satisfy non-tachyonic conditions at tree-level, skipping...')
                continue
            '''
            mcwd = bm.model(vh2, vs2, l1, l2, lm, ks2, yd, thetaY, m0, v2re)

            # Physical parameters
            phy = physpara(mcwd)
            if type(phy) is not list and phy == 1:
                ls_novh.append(l2)
                lsh_novh.append(lm)
                yd_novh.append(yd)
                m0_novh.append(m0)
                continue
            elif type(phy) is not list and phy == 2:
                ls_nomh.append(l2)
                lsh_nomh.append(lm)
                yd_nomh.append(yd)
                m0_nomh.append(m0)
            elif type(phy) is not list and phy == 3:
                ls_vs.append(l2)
                lsh_vs.append(lm)
                yd_vs.append(yd)
                m0_vs.append(m0)
                continue
            else:
            #vphy, tanbphy, m1phy, m2phy, sintphy = phy[0], phy[1], phy[2], phy[3], phy[4]
                vphy, m1phy, localMin = phy[0], phy[1], phy[2]
                ls_vh.append(l2)
                lsh_vh.append(lm)
                yd_vh.append(yd)
                m0_vh.append(m0)
                scan_rank.append([vh2, vs2, l1, l2, lm, ks2, yd, thetaY, m0, v2re, vphy, m1phy, localMin])
                filename = '%s_%s' % (FILE, rank)
                np.save(os.path.join(BASE_PATH, filename), scan_rank)

            #Tong: Constraints: 1. Higgs mass and vev.
            '''            
            if all((vphy <= 248., vphy >= 244.)):
                
                if all ((m1phy <= 127., m1phy >= 123.)):

                    scan_rank.append([m12, m22, l1, l2, lm, ks2, yd, thetaY, m0, v2re, vphy, m1phy, localMin])
                    filename = '%s_%s' % (FILE, rank)
                    np.save(os.path.join(BASE_PATH, filename), scan_rank)

                    
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
                        ks2_sfo.append(ks2)

                    if fo:
                        if not sfo:
                            yd_fo.append(yd)
                            # v/Tc
                            # sh_fo.append(foph[0])
                            ls_fo.append(l2)
                            lsh_fo.append(lm)
                            ks2_fo.append(ks2)
                    
                        scan_rank.append([m12, m22, l1, l2, lm, ks2, yd, thetaY, m0, v2re, vphy, tanbphy, m1phy, m2phy, sintphy])
                        filename = '%s_%s' % (FILE, rank)
                        np.save(os.path.join(BASE_PATH, filename), scan_rank)

                    else:
            
                        ls_scan.append(l2)
                        lsh_scan.append(lm)
                        yd_scan.append(yd)
                        ks2_scan.append(ks2)
                    
                else:
                                        
                    pass
           
            else:
                
                pass
    
    '''
    '''
    para_dict = {'ls_scan':ls_scan, 'lsh_scan':lsh_scan, 'yd_scan':yd_scan, 'ks2_scan':ks2_scan,
            'ls_fo':ls_fo, 'lsh_fo':lsh_fo, 'yd_fo':yd_fo, 'ks2_fo': ks2_fo,
            'ls_sfo':ls_sfo, 'lsh_sfo':lsh_sfo, 'yd_sfo':yd_sfo, 'ks2_sfo': ks2_sfo}
            '''
    para_dict = { 'ls_novh':ls_novh, 'lsh_novh':lsh_novh, 'ls_nomh':ls_nomh, 'lsh_nomh':lsh_nomh,
                'ls_vs':ls_vs, 'lsh_vs':lsh_vs, 'ls_vh':ls_vh, 'lsh_vh':lsh_vh,
                'yd_novh':yd_novh, 'm0_novh':m0_novh, 'yd_nomh':yd_nomh, 'm0_nomh':m0_nomh,
                'yd_vs':yd_vs, 'm0_vs':m0_vs, 'yd_vh':yd_vh, 'm0_vh':m0_vh}
    filename = '%s/%s_%s_plot.h5' % (BASE_PATH, FILE, rank)
    if os.path.exists(filename):
        os.remove(filename)
    dd.io.save(filename, para_dict)

                                                 
scan = getscan(lsr, lshr, ks2r, vsr, ydr, thetaYr, m0r, npt)
tf = time.clock()
dt = tf-t0
print ('Run time: %s' % dt)


