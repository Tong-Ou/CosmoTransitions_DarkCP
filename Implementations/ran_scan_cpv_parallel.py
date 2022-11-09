#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 17:10:18 2018

@author: Tong Ou (modified from Yikun's codes)
"""

"""
In this module, we create a scaning over the parameter space (ms, sint, tanb) which enters the model 
defined and used in the module baseMo.py

We create a dictionary containing the parameter info and the zero temperature local minima info. 

Usage functions are created to print/plot the scaning. 
"""

import sys
import baseMo_s_cpv as bm

import deepdish as dd

import random as ran

import numpy as np

import time

import os


#import mpi4py.MPI as mpi

import math

import matplotlib.pyplot as plt

import time

#comm = mpi.COMM_WORLD

rank = int(sys.argv[1]) #task_id

numtasks = int(sys.argv[2])

FILE = sys.argv[3]

npt = int(sys.argv[4])

t0 = time.clock()

#Fixed parameters
lh = 0.129
vh = 246.

BASE_PATH = os.getcwd()

#Bounds from T=0 boundary conditions
def ls_max(lh, vh2, vs2, ks2):
    if lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2 >= 0:
	ls_max = (lh*vh2**2-4.*ks2*vs2+np.sqrt(lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2))/(2.*vs2**2)
        return ls_max
    else:
	return lh*vh2**2/(2.*vs2**2)

def ls_min(lh, vh2, ma2, lsh):
    return (-2*ma2+lsh*vh2)**2./(4*lh*vh2**2)

#Check whether the higgs vacuum is the global minimum
def physpara(m):
    print ('\n')
    print ('ls:%s lsh:%s vs2:%s ks2:%s m0:%s yd:%s thetaYd:%s ' % (m.ls, m.lsh, m.vs2, m.ks2, m.m0, m.yd, m.thetaY))
    print ('non-tachyonic: %s' % bool((0.5*m.lsh*m.vh2 - 2*m.ks2)/m.vs2-m.ls>0.))
    print ('condition 2-a: %s' % all(((m.lh*m.vh2**2-4.*m.ks2*m.vs2+np.sqrt(m.lh**2*m.vh2**4-8.*m.lh*m.vh2**2*m.ks2*m.vs2))/(2.*m.vs2**2)-m.ls>0., (m.lh*m.vh2**2-4.*m.ks2*m.vs2-np.sqrt(m.lh**2*m.vh2**4-8.*m.lh*m.vh2**2*m.ks2*m.vs2))/(2.*m.vs2**2)-m.ls<0. )))

    vp = m.findMinimum(T=0.)
    print ('Higgs minimum: [%s, %s, %s]' % (vp[0],vp[1],vp[2]))
    if vp[0] > 248. or vp[0] < 244.:
        d2V_exp = m.d2V([246.,0.,0.], T=0.)
        m2sq_exp, eig_exp = np.linalg.eig(d2V_exp)
        print ('\n')
        print ('Mass eigenvalues at expected Higgs minimum: %s' % m2sq_exp)
	# Note that the eigenvalues returned by np.linalg.eig might not be sorted in a particular order!!!
        if any((m2sq_exp[0]<0, m2sq_exp[1]<0, m2sq_exp[2]<0)):
            print ('Expected Higgs vev is a saddle point, exitting...')
            return 1
	else:
            print ('Local minimum %s is deeper than expected Higgs minimum, exitting...' % vp)
            return 3
    else:
        m2phys = m.d2V(vp, T=0.)
        m2physd, eigv = np.linalg.eig(m2phys)
        # Minimum found by python is normally a real minimum, here is a double check.
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
            vs = abs((m.vs2+2.*m.ks2/m.ls))**0.5
            vscan = [0.1, vh, vs, 1e4]
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
   
            print ('T=0 boundary conditions met! Recording BM...')
            return [vp, m1phy, minX]
        
        
def getscan(npt):
       
    scan_task = range(npt)
    
    #rank = comm.Get_rank()

    #size = comm.Get_size()
        
    scan_rank, vphy_list, lmin_list = [], [], []
                
    ran.seed(time.time())

    ls_novh, lsh_novh, vs2_novh, ks2_novh, yd_novh, thetaY_novh, m0_novh = [], [], [], [], [], [], []
    ls_nomh, lsh_nomh, vs2_nomh, ks2_nomh, yd_nomh, thetaY_nomh, m0_nomh = [], [], [], [], [], [], []
    ls_vs, lsh_vs, vs2_vs, ks2_vs, yd_vs, thetaY_vs, m0_vs = [], [], [], [], [], [], []
    ls_vh, lsh_vh, vs2_vh, ks2_vh, yd_vh, thetaY_vh, m0_vh = [], [], [], [], [], [], []

    logfile = '%s/%s_%s.log' % (BASE_PATH, FILE, rank)
    log = open(logfile, 'w')
    sys.stdout = log 
        
    for n in range(len(scan_task[rank:npt:numtasks])):

	l1 = lh
	vh2 = vh**2.
	ma = ran.uniform(1, 10.)
	ma2 = ma**2.
	dm = ran.uniform(5., 35.)
	ms = ma + dm
        ks2 = (ms**2 - ma2)/4. #np.exp(ran.uniform(np.log(ks2min), np.log(ks2max)))
        thetaY = ran.uniform(0., np.pi/2)
	m0min = ms
        m0 = ran.uniform(m0min, 3.*m0min)
	yd = ran.uniform((m0/1.4e3)**0.5, (m0/1e3)**0.5)
	v2re = 172.9**2.

	vsp = ran.uniform(1.0e3, 2.5e3)

        lshmin = 2*ma2/vh2
        lshmax = 0.01 #min(0.01, 2*(l1*vh2/vsp**2 + ma2/vh2)*1.05)
        if lshmin < lshmax:
            lm = ran.uniform(lshmin, lshmax)
        else:
            continue

	l2 = (0.5*lm*vh2 - ma2)/(vsp**2)
	vs2 = vsp**2 - 2*ks2/l2
        
        mcwd = bm.model(vh2, vs2, l1, l2, lm, ks2, yd, thetaY, m0, v2re)

        # Physical parameters
        phy = physpara(mcwd)
        if type(phy) is not list and phy == 1:
            ls_novh.append(l2)
            lsh_novh.append(lm)
	    vs2_novh.append(vs2)
	    ks2_novh.append(ks2)
            yd_novh.append(yd)
	    thetaY_novh.append(thetaY)
            m0_novh.append(m0)
        elif type(phy) is not list and phy == 2:
            ls_nomh.append(l2)
            lsh_nomh.append(lm)
	    vs2_nomh.append(vs2)
	    ks2_nomh.append(ks2)
            yd_nomh.append(yd)
	    thetaY_nomh.append(thetaY)
            m0_nomh.append(m0)
        elif type(phy) is not list and phy == 3:
            ls_vs.append(l2)
            lsh_vs.append(lm)
	    vs2_vs.append(vs2)
	    ks2_vs.append(ks2)
            yd_vs.append(yd)
	    thetaY_vs.append(thetaY)
            m0_vs.append(m0)
        else:
        #vphy, tanbphy, m1phy, m2phy, sintphy = phy[0], phy[1], phy[2], phy[3], phy[4]
            vphy, m1phy, localMin = phy[0], phy[1], phy[2]
            ls_vh.append(l2)
            lsh_vh.append(lm)
	    vs2_vh.append(vs2)
	    ks2_vh.append(ks2)
            yd_vh.append(yd)
	    thetaY_vh.append(thetaY)
            m0_vh.append(m0)
            scan_rank.append([vh2, vs2, l1, l2, lm, ks2, yd, thetaY, m0, v2re, m1phy])
	    vphy_list.append(vphy)
	    lmin_list.append(localMin)
            filename = '%s_%s' % (FILE, rank)
	    vphy_file = filename + '_vphy'
	    lmin_file = filename + '_localMin'
            np.save(os.path.join(BASE_PATH, filename), scan_rank, allow_pickle=True)
	    np.save(os.path.join(BASE_PATH, vphy_file), vphy_list, allow_pickle=True)
	    np.save(os.path.join(BASE_PATH, lmin_file), lmin_list, allow_pickle=True)

        para_dict = { 'ls_novh':ls_novh, 'lsh_novh':lsh_novh, 'ls_nomh':ls_nomh, 'lsh_nomh':lsh_nomh,
                        'ls_vs':ls_vs, 'lsh_vs':lsh_vs, 'ls_vh':ls_vh, 'lsh_vh':lsh_vh,
			'vs2_novh':vs2_novh, 'ks2_novh':ks2_novh, 'vs2_nomh':vs2_nomh, 'ks2_nomh':ks2_nomh,
			'vs2_vs':vs2_vs, 'ks2_vs':ks2_vs, 'vs2_vh':vs2_vh, 'ks2_vh':ks2_vh,
                        'yd_novh':yd_novh, 'thetaY_novh':thetaY_novh, 'm0_novh':m0_novh, 
			'yd_nomh':yd_nomh, 'thetaY_nomh':thetaY_nomh, 'm0_nomh':m0_nomh,
                        'yd_vs':yd_vs, 'thetaY_vs':thetaY_vs, 'm0_vs':m0_vs,
			 'yd_vh':yd_vh, 'thetaY_vh':thetaY_vh, 'm0_vh':m0_vh}
        filename = '%s/%s_%s_plot.h5' % (BASE_PATH, FILE, rank)
        dd.io.save(filename, para_dict)

if __name__ == "__main__":   
    getscan(npt)
    tf = time.clock()
    dt = tf-t0
    print ('Run time: %s' % dt)


