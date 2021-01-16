'''
author: Tong Ou
This script is to plot the viable parameter space based on results of ran_scan.py.
'''

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import sys 
import deepdish as dd
import argparse
import os

#ncpu = 4

ls_scan = np.zeros(0)
ls_sfo = np.zeros(0)
lsh_scan = np.zeros(0)
lsh_sfo = np.zeros(0)
sh_fo = np.zeros(0)
yd_fo = np.zeros(0)

parser = argparse.ArgumentParser()
parser.add_argument('--FILE', type=str, default='', help='File name prefix.')
parser.add_argument('--fullScan', dest='fullScan', default=False, action='store_true', help='Whether scanning over the entire parameter space or not.')

args = parser.parse_args()
FILE = args.FILE
fullScan = args.fullScan

def ls_max_nontach(lsh, vh2, vs2, ks2):
    return (0.5*lsh*vh2 - 2*ks2)/vs2

def ls_max_IIc(lsh, lh):
    return lsh**2./(4*lh)

def ls_range_IIa(lh, vh2, vs2, ks2):
    if lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2 >=0:
        ls_min = (lh*vh2**2-4.*ks2*vs2-np.sqrt(lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2))/(2.*vs2**2)
        ls_max = (lh*vh2**2-4.*ks2*vs2+np.sqrt(lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2))/(2.*vs2**2)
        return [ls_min, ls_max]
    else:
        return None

def plot_min(FILE, ncpu, lh, vh2, vs2, ks2, xlim, ylim):

    '''
    for i in range(ncpu):
        para_dict = dd.io.load('%s_%s_plot.h5' % (FILE, i))
        ls_novh_i, ls_nomh_i, ls_vs_i, ls_vh_i = para_dict['ls_novh'], para_dict['ls_nomh'], para_dict['ls_vs'], para_dict['ls_vh']
        lsh_nomh_i, lsh_novh_i, lsh_vs_i, lsh_vh_i = para_dict['lsh_nomh'], para_dict['lsh_novh'], para_dict['lsh_vs'], para_dict['lsh_vh']
        yd_novh_i, yd_nomh_i, yd_vs_i, yd_vh_i = para_dict['yd_novh'], para_dict['yd_nomh'], para_dict['yd_vs'], para_dict['yd_vh']
        m0_novh_i, m0_nomh_i, m0_vs_i, m0_vh_i = para_dict['m0_novh'], para_dict['m0_nomh'], para_dict['m0_vs'], para_dict['m0_vh']

        if i == 0:
            ls_nomh, ls_novh, ls_vs, ls_vh = ls_nomh_i, ls_novh_i, ls_vs_i, ls_vh_i
            lsh_nomh, lsh_novh, lsh_vs, lsh_vh = lsh_nomh_i, lsh_novh_i, lsh_vs_i, lsh_vh_i
            yd_novh, yd_nomh, yd_vs, yd_vh = yd_novh_i, yd_nomh_i, yd_vs_i, yd_vh_i
            m0_novh, m0_nomh, m0_vs, m0_vh = m0_novh_i, m0_nomh_i, m0_vs_i, m0_vh_i
        else:
            ls_nomh = np.concatenate((ls_nomh, ls_nomh_i), axis=0)
            ls_novh = np.concatenate((ls_novh, ls_novh_i), axis=0)
            ls_vs = np.concatenate((ls_vs, ls_vs_i), axis=0)
            ls_vh = np.concatenate((ls_vh, ls_vh_i), axis=0)
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
	'''
    para_dict = dd.io.load('%s_plot.h5' % FILE)
    ls_novh, ls_nomh, ls_vs, ls_vh = para_dict['ls_novh'], para_dict['ls_nomh'], para_dict['ls_vs'], para_dict['ls_vh']
    lsh_nomh, lsh_novh, lsh_vs, lsh_vh = para_dict['lsh_nomh'], para_dict['lsh_novh'], para_dict['lsh_vs'], para_dict['lsh_vh']
    yd_novh, yd_nomh, yd_vs, yd_vh = para_dict['yd_novh'], para_dict['yd_nomh'], para_dict['yd_vs'], para_dict['yd_vh']
    m0_novh, m0_nomh, m0_vs, m0_vh = para_dict['m0_novh'], para_dict['m0_nomh'], para_dict['m0_vs'], para_dict['m0_vh']    

    # Save the plotting variables to a single file
    '''
    para_dict = { 'ls_novh':ls_novh, 'lsh_novh':lsh_novh, 'ls_nomh':ls_nomh, 'lsh_nomh':lsh_nomh,
                'ls_vs':ls_vs, 'lsh_vs':lsh_vs, 'ls_vh':ls_vh, 'lsh_vh':lsh_vh,
                'yd_novh':yd_novh, 'm0_novh':m0_novh, 'yd_nomh':yd_nomh, 'm0_nomh':m0_nomh,
                'yd_vs':yd_vs, 'm0_vs':m0_vs, 'yd_vh':yd_vh, 'm0_vh':m0_vh}
    filename = '%s_plot.h5' % FILE
    if os.path.exists(filename):
        os.remove(filename)
    dd.io.save(filename, para_dict)
    '''

    ls_min, ls_max = None, None
    ls_lim = ls_range_IIa(lh, vh2, vs2, ks2)
    if ls_lim is not None:
        ls_min = ls_lim[0] * np.ones(100)
        ls_max = ls_lim[1] * np.ones(100)

    fig = plt.figure(figsize=(8,6))
    ax1 = plt.subplot(111)
    lshr = np.linspace(0., 4*np.pi, 100)
    ax1.plot(lshr, ls_max_nontach(lshr, vh2, vs2, ks2), c='black')
    #ax1.plot(lshr, ls_max_IIc(lshr, lh), c='magenta')
    if ls_min is not None:
        ax1.plot(lshr, ls_min, c='blue')
        ax1.plot(lshr, ls_max, c='blue')
    size = 10.
    ax1.scatter(lsh_novh, ls_novh, s=size, c='black', alpha=0.5, label='false Higgs vev')
    ax1.scatter(lsh_nomh, ls_nomh, s=size, c='green', alpha=0.5, label='false Higgs mass')
    ax1.scatter(lsh_vs, ls_vs, s=size, c='blue', alpha=0.5, label='false minimum')
    ax1.scatter(lsh_vh, ls_vh, s=size, c='red', label='meet BC')
    ax1.legend(fontsize=15, loc='upper right')
    ax1.set_ylim(0., ylim)
    ax1.set_xlim(0., xlim)
    ax1.set_xlabel(r'$\lambda_{SH}$', size=15)
    ax1.set_ylabel(r'$\lambda_S$', size=15)
    plt.savefig('%s_lsh_ls.pdf' % (FILE))

    fig2 = plt.figure(figsize=(8,6))
    ax2 = plt.subplot(111)
    ax2.scatter(m0_novh, yd_novh, s=size, c='black', alpha=0.5, label='false Higgs vev')
    ax2.scatter(m0_nomh, yd_nomh, s=size, c='green', alpha=0.5, label='false Higgs mass')
    ax2.scatter(m0_vs, yd_vs, s=size, c='blue', alpha=0.5, label='false minimum')
    ax2.scatter(m0_vh, yd_vh, s=size, c='red', label='meet BC')
    ax2.legend(fontsize=15, loc='upper right')
    ax2.set_xlabel(r'$m_0\ [GeV]$', size=15)
    ax2.set_ylabel(r'$\lambda$', size=15)
    plt.savefig('%s_m0_yd.pdf' % (FILE))

def plot_min_ana(lh, vh2, vs2, ks2, xlim, ylim):

    ls_min, ls_max = None, None
    ls_lim = ls_range_IIa(lh, vh2, vs2, ks2)
    if ls_lim is not None:
        ls_min = ls_lim[0] * np.ones(100)
        ls_max = ls_lim[1] * np.ones(100)

    fig = plt.figure(figsize=(8,6))
    ax1 = plt.subplot(111)
    lshr = np.linspace(0., 4*np.pi, 100)
    ax1.plot(lshr, ls_max_nontach(lshr, vh2, vs2, ks2), c='black')
    #ax1.plot(lshr, ls_max_IIc(lshr, lh), c='magenta')
    ax1.fill_between(lshr, ls_max_nontach(lshr, vh2, vs2, ks2), ylim, facecolor='black', alpha=0.5, label='fail 1')
    #ax1.fill_between(lshr, ls_max_IIc(lshr, lh), ylim, facecolor='magenta', alpha=0.5, label='fail II-c')
    if ls_min is not None:
        ax1.plot(lshr, ls_min, c='blue')
        ax1.plot(lshr, ls_max, c='blue')
        ax1.fill_between(lshr, ls_max, ylim, facecolor='blue', alpha=0.5, label='fail 2-a')
        ax1.fill_between(lshr, 0., ls_min, facecolor='blue', alpha=0.5)
    else:
        ax1.fill_between(lshr, 0., ylim, facecolor='blue', alpha=0.5, label='fail 2-a')
    ax1.legend(fontsize=15)
    ax1.set_ylim(0., ylim)
    ax1.set_xlim(0., xlim)
    ax1.set_xlabel(r'$\lambda_{SH}$', size=15)
    ax1.set_ylabel(r'$\lambda_S$', size=15)
    plt.savefig('%s_lsh_ls_ana.pdf' % (FILE))

def plot_sfo(FILE, ncpu):
    for i in range(ncpu):
        para_dict = dd.io.load('%s_%s_plot.h5' % (FILE, i))
        ls_scan_i, ls_fo_i, ls_sfo_i = para_dict['ls_scan'], para_dict['ls_fo'], para_dict['ls_sfo']
        lsh_scan_i, lsh_fo_i, lsh_sfo_i = para_dict['lsh_scan'], para_dict['lsh_fo'], para_dict['lsh_sfo']
        yd_scan_i, yd_fo_i, yd_sfo_i = para_dict['yd_scan'], para_dict['yd_fo'], para_dict['yd_sfo']
    
        if i == 0:
            ls_scan, ls_fo, ls_sfo = ls_scan_i, ls_fo_i, ls_sfo_i
            lsh_scan, lsh_fo, lsh_sfo = lsh_scan_i, lsh_fo_i, lsh_sfo_i
            yd_scan, yd_fo, yd_sfo = yd_scan_i, yd_fo_i, yd_sfo_i
        else:
            ls_scan = np.concatenate((ls_scan, ls_scan_i), axis=0)
            ls_fo = np.concatenate((ls_fo, ls_fo_i), axis=0)
            ls_sfo = np.concatenate((ls_sfo, ls_sfo_i), axis=0)
            lsh_scan = np.concatenate((lsh_scan, lsh_scan_i), axis=0)
            lsh_fo = np.concatenate((lsh_fo, lsh_fo_i), axis=0)
            lsh_sfo = np.concatenate((lsh_sfo, lsh_sfo_i), axis=0)
            yd_scan = np.concatenate((yd_scan, yd_scan_i), axis=0)
            yd_fo = np.concatenate((yd_fo, yd_fo_i), axis=0)
            yd_sfo = np.concatenate((yd_sfo, yd_sfo_i), axis=0)
    
    lh = 0.129
    vh2 = 246.**2 #pay attention to the dtype!!
    vs2 = 150.**2
    ls = 2.
    ch = (9*0.65**2+3*0.35**2+2*(6*0.9911**2+12*lh))/48
    
    def yd_max(lsh):
        return 12*((ls/lh)*(vs2/vh2)*(ch+lsh/12.+lh)-(lh*ls)**0.5)-2*lsh-4*ls
    
    plt.figure(figsize=(12,5))
    ax = plt.subplot(121)
    ax.scatter(ls_scan, lsh_scan, c='black', label='Scan Points')
    ax.scatter(ls_fo, lsh_fo, c='blue', label='First-order Trans')
    ax.scatter(ls_sfo, lsh_sfo, c='red', label=r'$v_c/T_c>1$')
    ax.legend()
    ax.set_xlabel(r'$\lambda_S$')
    ax.set_ylabel(r'$\lambda_{SH}$')
    ax.grid(True)
    
    ax = plt.subplot(122)
    #ax.scatter(yd_fo, sh_fo, c='blue', label='First-order Trans')
    #ax.scatter(yd_sfo, sh_sfo, c='red', label='sfo')
    #ax.legend()
    #ax.set_xlabel(r'$y_d$')
    #ax.set_ylabel(r'$v_c/T_c$')
    #ax.grid(True)
    #
    #plt.show()    
    
    #plt.figure()
    lsh = np.linspace(0, 4*np.pi, 100)
    plt.scatter(lsh_scan, np.square(yd_scan), c='black', label='Scan Points')
    plt.scatter(lsh_fo, np.square(yd_fo), c='blue', label='First-order Trans')
    plt.scatter(lsh_sfo, np.square(yd_sfo), c='red', label=r'$v_c/T_c>1$')
    if not fullScan:
        plt.plot(lsh, yd_max(lsh), linewidth=2., c='black', label=r'Upper bound for $|y_d|^2$ from high-T')
    #plt.ylim(0,5)
    plt.legend()
    plt.xlabel(r'$\lambda_{SH}$')
    plt.ylabel(r'$|y_d|^2$')
    plt.show()

#===================================================
#Main function
#===================================================
lh = 0.129
vh2 = 246.**2
vs2 = 100.**2
ks2 = 3000.

plot_min(FILE, 4, lh, vh2, vs2, ks2, 4*np.pi, 6.)
plot_min_ana(lh, vh2, vs2, ks2, 4*np.pi, 6.)

