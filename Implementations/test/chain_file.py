'''
This script is to chain the files produced by multiple 
CPU cores.
'''
import numpy as np
import os
import sys
import deepdish as dd

ncpu = int(sys.argv[1])
FILE = sys.argv[2]
#FILE2 = sys.argv[2]
#FILE3 = sys.argv[3]
para_list = []
vphy_list = []
lmin_list = []

for i in range(ncpu):
#for filename in [FILE, FILE2]:
    filename = '%s_%s' % (FILE, i)
    try:
        para = np.load(filename + '.npy', allow_pickle = True).tolist()
        para_list.append(para)
	vphy = np.load(filename + '_vphy.npy', allow_pickle = True).tolist()
	vphy_list.append(vphy)
	lmin = np.load(filename + '_localMin.npy', allow_pickle = True).tolist()
	lmin_list.append(lmin)
    except:
        continue
 
if len(para_list) > 0:
    paras = np.concatenate(para_list,axis = 0)
    np.save('%s.npy' % FILE, paras)
    vphys = np.concatenate(vphy_list,axis = 0)
    np.save('%s_vphy.npy' % FILE, vphys)
    lmins = np.concatenate(lmin_list,axis = 0)
    np.save('%s_localMin.npy' % FILE, lmins)
else:
    pass

for i in range(ncpu):
    try:
    	para_dict = dd.io.load('%s_%s_plot.h5' % (FILE, i))
    except:
	continue
    ls_novh_i, ls_nomh_i, ls_vs_i, ls_vh_i = para_dict['ls_novh'], para_dict['ls_nomh'], para_dict['ls_vs'], para_dict['ls_vh']
    lsh_nomh_i, lsh_novh_i, lsh_vs_i, lsh_vh_i = para_dict['lsh_nomh'], para_dict['lsh_novh'], para_dict['lsh_vs'], para_dict['lsh_vh']
    vs2_novh_i, vs2_nomh_i, vs2_vs_i, vs2_vh_i = para_dict['vs2_novh'], para_dict['vs2_nomh'], para_dict['vs2_vs'], para_dict['vs2_vh']
    ks2_novh_i, ks2_nomh_i, ks2_vs_i, ks2_vh_i = para_dict['ks2_novh'], para_dict['ks2_nomh'], para_dict['ks2_vs'], para_dict['ks2_vh']
    yd_novh_i, yd_nomh_i, yd_vs_i, yd_vh_i = para_dict['yd_novh'], para_dict['yd_nomh'], para_dict['yd_vs'], para_dict['yd_vh']
    m0_novh_i, m0_nomh_i, m0_vs_i, m0_vh_i = para_dict['m0_novh'], para_dict['m0_nomh'], para_dict['m0_vs'], para_dict['m0_vh']
    thetaY_novh_i, thetaY_nomh_i, thetaY_vs_i, thetaY_vh_i = para_dict['thetaY_novh'], para_dict['thetaY_nomh'], para_dict['thetaY_vs'], para_dict['thetaY_vh']

    if i == 0:
        ls_nomh, ls_novh, ls_vs, ls_vh = ls_nomh_i, ls_novh_i, ls_vs_i, ls_vh_i
        lsh_nomh, lsh_novh, lsh_vs, lsh_vh = lsh_nomh_i, lsh_novh_i, lsh_vs_i, lsh_vh_i
	vs2_nomh, vs2_novh, vs2_vs, vs2_vh = vs2_nomh_i, vs2_novh_i, vs2_vs_i, vs2_vh_i
	ks2_nomh, ks2_novh, ks2_vs, ks2_vh = ks2_nomh_i, ks2_novh_i, ks2_vs_i, ks2_vh_i
        yd_novh, yd_nomh, yd_vs, yd_vh = yd_novh_i, yd_nomh_i, yd_vs_i, yd_vh_i
        m0_novh, m0_nomh, m0_vs, m0_vh = m0_novh_i, m0_nomh_i, m0_vs_i, m0_vh_i
	thetaY_novh, thetaY_nomh, thetaY_vs, thetaY_vh = thetaY_novh_i, thetaY_nomh_i, thetaY_vs_i, thetaY_vh_i
    else:
        ls_nomh = np.concatenate((ls_nomh, ls_nomh_i), axis=0)
        ls_novh = np.concatenate((ls_novh, ls_novh_i), axis=0)
        ls_vs = np.concatenate((ls_vs, ls_vs_i), axis=0)
        ls_vh = np.concatenate((ls_vh, ls_vh_i), axis=0)

        vs2_nomh = np.concatenate((vs2_nomh, vs2_nomh_i), axis=0)
        vs2_novh = np.concatenate((vs2_novh, vs2_novh_i), axis=0)
        vs2_vs = np.concatenate((vs2_vs, vs2_vs_i), axis=0)
        vs2_vh = np.concatenate((vs2_vh, vs2_vh_i), axis=0)

        ks2_nomh = np.concatenate((ks2_nomh, ks2_nomh_i), axis=0)
        ks2_novh = np.concatenate((ks2_novh, ks2_novh_i), axis=0)
        ks2_vs = np.concatenate((ks2_vs, ks2_vs_i), axis=0)
        ks2_vh = np.concatenate((ks2_vh, ks2_vh_i), axis=0)

        lsh_nomh = np.concatenate((lsh_nomh, lsh_nomh_i), axis=0)
        lsh_novh = np.concatenate((lsh_novh, lsh_novh_i), axis=0)
        lsh_vs = np.concatenate((lsh_vs, lsh_vs_i), axis=0)
        lsh_vh = np.concatenate((lsh_vh, lsh_vh_i), axis=0)

        yd_nomh = np.concatenate((yd_nomh, yd_nomh_i), axis=0)
        yd_novh = np.concatenate((yd_novh, yd_novh_i), axis=0)
        yd_vs = np.concatenate((yd_vs, yd_vs_i), axis=0)
        yd_vh = np.concatenate((yd_vh, yd_vh_i), axis=0)

        m0_nomh = np.concatenate((m0_nomh, m0_nomh_i), axis=0)
        m0_novh = np.concatenate((m0_novh, m0_novh_i), axis=0)
        m0_vs = np.concatenate((m0_vs, m0_vs_i), axis=0)
        m0_vh = np.concatenate((m0_vh, m0_vh_i), axis=0)

        thetaY_nomh = np.concatenate((thetaY_nomh, thetaY_nomh_i), axis=0)
        thetaY_novh = np.concatenate((thetaY_novh, thetaY_novh_i), axis=0)
        thetaY_vs = np.concatenate((thetaY_vs, thetaY_vs_i), axis=0)
        thetaY_vh = np.concatenate((thetaY_vh, thetaY_vh_i), axis=0)

para_dict = { 'ls_novh':ls_novh, 'lsh_novh':lsh_novh, 'ls_nomh':ls_nomh, 'lsh_nomh':lsh_nomh,
            'ls_vs':ls_vs, 'lsh_vs':lsh_vs, 'ls_vh':ls_vh, 'lsh_vh':lsh_vh,
	    'vs2_novh':vs2_novh, 'ks2_novh':ks2_novh, 'vs2_nomh':vs2_nomh, 'ks2_nomh':ks2_nomh,
	    'vs2_vs':vs2_vs, 'ks2_vs':ks2_vs, 'vs2_vh':vs2_vh, 'ks2_vh':ks2_vh,
            'yd_novh':yd_novh, 'thetaY_novh':thetaY_novh, 'm0_novh':m0_novh, 'yd_nomh':yd_nomh, 'thetaY_nomh':thetaY_nomh, 'm0_nomh':m0_nomh,
            'yd_vs':yd_vs, 'm0_vs':m0_vs, 'thetaY_vs':thetaY_vs, 'yd_vh':yd_vh, 'm0_vh':m0_vh, 'thetaY_vh':thetaY_vh}
filename = '%s_plot.h5' % FILE
if os.path.exists(filename):
    os.remove(filename)
dd.io.save(filename, para_dict)

