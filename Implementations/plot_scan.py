'''
author: Tong Ou
This script is to plot the viable parameter space based on results of ran_scan.py.
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

    # Load data
    para_dict = dd.io.load('%s_plot.h5' % FILE)
    ls_novh, ls_nomh, ls_vh = para_dict['ls_novh'], para_dict['ls_nomh'], para_dict['ls_vh']
    lsh_nomh, lsh_novh, lsh_vh = para_dict['lsh_nomh'], para_dict['lsh_novh'], para_dict['lsh_vh']
    yd_novh, yd_nomh, yd_vh = para_dict['yd_novh'], para_dict['yd_nomh'], para_dict['yd_vh']
    m0_novh, m0_nomh, m0_vh = para_dict['m0_novh'], para_dict['m0_nomh'], para_dict['m0_vh']    

    para_vs = dd.io.load('%s_vs.h5' % FILE)
    ls_inf, ls_vs, ls_va, ls_sa, ls_as, ls_orig, ls_none = para_vs['ls_inf'], para_vs['ls_vs'], para_vs['ls_va'], para_vs['ls_sa'], para_vs['ls_as'], para_vs['ls_orig'], para_vs['ls_none']
    lsh_inf, lsh_vs, lsh_va, lsh_sa, lsh_as, lsh_orig, lsh_none = para_vs['lsh_inf'], para_vs['lsh_vs'], para_vs['lsh_va'], para_vs['lsh_sa'], para_vs['lsh_as'], para_vs['lsh_orig'], para_vs['lsh_none']
    m0_inf, m0_vs, m0_va, m0_sa, m0_as, m0_orig, m0_none = para_vs['m0_inf'], para_vs['m0_vs'], para_vs['m0_va'], para_vs['m0_sa'], para_vs['m0_as'], para_vs['m0_orig'], para_vs['m0_none']
    yd_inf, yd_vs, yd_va, yd_sa, yd_as, yd_orig, yd_none = para_vs['yd_inf'], para_vs['yd_vs'], para_vs['yd_va'], para_vs['yd_sa'], para_vs['yd_as'], para_vs['yd_orig'], para_vs['yd_none']

    print ('Total scan rounds:')
    print (len(ls_novh)+len(ls_nomh)+len(ls_vh)+len(ls_inf)+len(ls_vs)+len(ls_va)+len(ls_sa)+len(ls_as)+len(ls_orig)+len(ls_none))
    # Draw tree-level boundaries
    ls_min, ls_max = None, None
    ls_lim = ls_range_IIa(lh, vh2, vs2, ks2)
    if ls_lim is not None:
        ls_min = ls_lim[0] * np.ones(100)
        ls_max = ls_lim[1] * np.ones(100)

    fig = plt.figure(figsize=(12,6))
    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2,1])
    ax1 = fig.add_subplot(spec[0,0])
    lshr = np.linspace(0., 4*np.pi, 100)
    
    
    ax1.plot(lshr, ls_max_nontach(lshr, vh2, vs2, ks2), c='black')
    #ax1.plot(lshr, ls_max_IIc(lshr, lh), c='magenta')
    if ls_min is not None:
        ax1.plot(lshr, ls_min, c='blue')
        ax1.plot(lshr, ls_max, c='blue')
    
    # Plot args
    size = 10.
    evenly_spaced_interval = range(7)
    colors = [plt.cm.Dark2(x) for x in evenly_spaced_interval]
    
    ax1.scatter(lsh_novh, ls_novh, s=size, c='grey', alpha=0.5, label='false Higgs vev')
    ax1.scatter(lsh_nomh, ls_nomh, s=size, c='black', alpha=0.8, label='false Higgs mass')
    ax1.scatter(lsh_inf, ls_inf, s=size, c=colors[0], alpha=0.5, label='Potential not bounded from below')
    ax1.scatter(lsh_va, ls_va, s=size, c=colors[1], alpha=0.8, label=r'$V(0,0,a)<V(h,0,0), |a|>0$')
    ax1.scatter(lsh_vs, ls_vs, s=size, c=colors[2], alpha=0.8, label=r'$V(0,s,0)<V(h,0,0), |s|>0$')
    ax1.scatter(lsh_sa, ls_sa, s=size, c=colors[3], alpha=0.8, label=r'$V(0,s,a)<V(h,0,0), |s|>|a|>0$')
    ax1.scatter(lsh_as, ls_as, s=size, c=colors[4], alpha=0.5, label=r'$V(0,s,a)<V(h,0,0), |a|>|s|>0$')
    ax1.scatter(lsh_orig, ls_orig, s=size, c=colors[5], alpha=0.8, label=r'$V(0,0,0)<V(h,0,0)$')
    ax1.scatter(lsh_none, ls_none, s=size, c=colors[6], alpha=0.8, label=r'false global minimum')
    ax1.scatter(lsh_vh, ls_vh, s=size, c='red', label='meet BC')
    ax1.legend(fontsize=10, loc='upper right', bbox_to_anchor=(1.55, 1))
    ax1.set_ylim(0., ylim)
    ax1.set_xlim(0., xlim)
    ax1.set_xlabel(r'$\lambda_{SH}$', size=15)
    ax1.set_ylabel(r'$\lambda_S$', size=15)
    plt.savefig('%s_lsh_ls_2.pdf' % (FILE))

    
    fig2 = plt.figure(figsize=(12,6))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2,1])
    ax2 = fig2.add_subplot(spec2[0,0])
    ax2.scatter(m0_novh, yd_novh, s=size, c='black', alpha=0.5, label='false Higgs vev')
    ax2.scatter(m0_nomh, yd_nomh, s=size, c='grey', alpha=0.5, label='false Higgs mass')
    ax2.scatter(m0_inf, yd_inf, s=size, c=colors[0], alpha=0.5, label='Potential not bounded from below')
    ax2.scatter(m0_va, yd_va, s=size, c=colors[1], alpha=0.5, label=r'$V(0,0,a)<V(h,0,0), |a|>0$')
    ax2.scatter(m0_vs, yd_vs, s=size, c=colors[2], alpha=0.5, label=r'$V(0,s,0)<V(h,0,0), |s|>0$')
    ax2.scatter(m0_sa, yd_sa, s=size, c=colors[3], alpha=0.7, label=r'$V(0,s,a)<V(h,0,0), |s|>|a|>0$')
    ax2.scatter(m0_as, yd_as, s=size, c=colors[4], alpha=0.7, label=r'$V(0,s,a)<V(h,0,0), |a|>|s|>0$')
    ax2.scatter(m0_orig, yd_orig, s=size, c=colors[5], alpha=0.5, label=r'$V(0,0,0)<V(h,0,0)$')
    ax2.scatter(m0_none, yd_none, s=size, c=colors[6], alpha=0.5, label=r'false global minimum')
    ax2.scatter(m0_vh, yd_vh, s=size, c='red', label='meet BC')
    ax2.legend(fontsize=10, loc='upper right', bbox_to_anchor=(1.55, 1))
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

#===================================================
#Main function
#===================================================
lh = 0.129
vh2 = 246.**2
vs2 = 100.**2
ks2 = 3000.

plot_min(FILE, 4, lh, vh2, vs2, ks2, 4*np.pi, 7.)
#plot_min_ana(lh, vh2, vs2, ks2, 4*np.pi, 4*np.pi)
