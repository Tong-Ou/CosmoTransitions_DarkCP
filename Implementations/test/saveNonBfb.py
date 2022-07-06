'''
author: Tong Ou
This script is to plot the phase transition calculation results.
'''

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys 
import deepdish as dd
import argparse
import os
import baseMo_s_cpv as bmt


rank = int(sys.argv[1]) #task_id
numtasks = int(sys.argv[2])
FILE = sys.argv[3]
filename = '%s_%s' % (FILE, rank)
para_dict = dd.io.load('%s_vs.h5' % FILE)

ls_inf, lsh_inf, vs2_inf, ks2_inf, m0_inf, yd_inf, thetaY_inf = para_dict['ls_inf'], para_dict['lsh_inf'], para_dict['vs2_inf'], para_dict['ks2_inf'], para_dict['m0_inf'], para_dict['yd_inf'], para_dict['thetaY_inf']

paras = []
vphys = []
lmins = []
scan_task = range(len(ls_inf))
rank_task = scan_task[rank:len(scan_task):numtasks]

for i in rank_task:
    vh2 = 246.**2
    lh = 0.129
    v2re = 172.9**2
    para = [vh2, vs2_inf[i], lh, ls_inf[i], lsh_inf[i], ks2_inf[i], yd_inf[i], thetaY_inf[i], m0_inf[i], v2re, 125.]
    m = bmt.model(vh2, vs2_inf[i], lh, ls_inf[i], lsh_inf[i], ks2_inf[i], yd_inf[i], thetaY_inf[i], m0_inf[i], v2re)
    vsp = abs((m.vs2+2.*m.ks2/m.ls))**0.5 
    vp = m.findMinimum(np.array([246., 0., 0.]), T=0.)
    if vp[0] > 248. or vp[0] < 244.:
	continue 
    vscan = [1e-5, 246., vsp]
    minXcand = []
    minX = []
    for v1 in vscan:
        for v2 in vscan:
            for v3 in vscan:
		minXcand.append([v1, v2, v3])

    vhGlobalMin = True
    for point in minXcand:
	lmin = m.findMinimum(point, T=0.)
	if all((abs(lmin[0]) < 5e5, abs(lmin[1]) < 5e5, abs(lmin[2]) < 5e5)):
	    minX.append(lmin)
	    if m.Vtot(lmin, T=0.) < m.Vtot(vp, T=0.):
	        if np.sum((abs(lmin)-abs(vp))**2 ,-1)**.5 >= 0.01:
		    print ('Local minimum %s is lower than Higgs minimum' % lmin)
		    print ('Higgs minimum is not global minimum, exiting...')
		    vhGlobalMin = False
		    break

    if not vhGlobalMin:
        continue                       
    else:
	paras.append(para)
	vphys.append(vp)
	lmins.append(minX)
	
	np.save(filename, paras)
	np.save(filename + '_vphy', vphys)
	np.save(filename + '_localMin', lmins)


