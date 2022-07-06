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

FILE_1 = sys.argv[1]
FILE_2 = sys.argv[2]
OUTDIR = FILE_1.rsplit('/', 1)[0]

def ls_max_nontach(lsh, vh2, vs2, ks2):
    return (0.5*lsh*vh2 - 2*ks2)/vs2

def ls_max_IIc(lsh, lh):
    return lsh**2./(4*lh)
'''
def ls_range_IIa(lh, vh2, vs2, ks2):
    if lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2 >=0:
        ls_min = (lh*vh2**2-4.*ks2*vs2-np.sqrt(lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2))/(2.*vs2**2)
        ls_max = (lh*vh2**2-4.*ks2*vs2+np.sqrt(lh**2*vh2**4-8.*lh*vh2**2*ks2*vs2))/(2.*vs2**2)
        return [ls_min, ls_max]
    else:
        return None
'''
def ls_min_IIa(lsh):
    ma2 = 10.**2
    lh = 0.129
    vh2 = 246.**2
    return (-2*ma2 + lsh*vh2)**2./(4*lh*vh2**2)

def lsh_min_bfb(vsp, m0, yd, ma):
    vh2 = 246.**2
    mt = 172.9
    md2 = m0**2+yd**2*10**8+np.sqrt(2)*m0*yd*10**4
    return 2/vh2*(yd**4*vsp**2/(16*np.pi**2)*(np.log(md2/mt**2)-3/2)+ma**2)

def lsh_max_IIa(vsp, ma):
    vh2 = 246.**2
    lh = 0.129
    return 2*(lh*vh2/vsp**2 + ma**2/vh2)

def plot_min(FILE, ma, ax, xmin, xmax, ymax):

    m0, yd = 80., 0.25
    
    # Load data
    para_dict = dd.io.load('%s_plot.h5' % FILE)
    ls_novh, ls_nomh, ls_vh = para_dict['ls_novh'], para_dict['ls_nomh'], para_dict['ls_vh']
    lsh_nomh, lsh_novh, lsh_vh = para_dict['lsh_nomh'], para_dict['lsh_novh'], para_dict['lsh_vh']
    vs2_novh, vs2_nomh, vs2_vh = para_dict['vs2_novh'], para_dict['vs2_nomh'], para_dict['vs2_vh']
    ks2_novh, ks2_nomh, ks2_vh = para_dict['ks2_novh'], para_dict['ks2_nomh'], para_dict['ks2_vh']
    yd_novh, yd_nomh, yd_vh = para_dict['yd_novh'], para_dict['yd_nomh'], para_dict['yd_vh']
    m0_novh, m0_nomh, m0_vh = para_dict['m0_novh'], para_dict['m0_nomh'], para_dict['m0_vh']    

    para_vs = dd.io.load('%s_vs.h5' % FILE)
    ls_inf, ls_vs, ls_va, ls_sa, ls_as, ls_orig, ls_none = para_vs['ls_inf'], para_vs['ls_vs'], para_vs['ls_va'], para_vs['ls_sa'], para_vs['ls_as'], para_vs['ls_orig'], para_vs['ls_none']
    lsh_inf, lsh_vs, lsh_va, lsh_sa, lsh_as, lsh_orig, lsh_none = para_vs['lsh_inf'], para_vs['lsh_vs'], para_vs['lsh_va'], para_vs['lsh_sa'], para_vs['lsh_as'], para_vs['lsh_orig'], para_vs['lsh_none']
    vs2_inf, vs2_vs, vs2_va, vs2_sa, vs2_as, vs2_orig, vs2_none = para_vs['vs2_inf'], para_vs['vs2_vs'], para_vs['vs2_va'], para_vs['vs2_sa'], para_vs['vs2_as'], para_vs['vs2_orig'], para_vs['vs2_none']
    ks2_inf, ks2_vs, ks2_va, ks2_sa, ks2_as, ks2_orig, ks2_none = para_vs['ks2_inf'], para_vs['ks2_vs'], para_vs['ks2_va'], para_vs['ks2_sa'], para_vs['ks2_as'], para_vs['ks2_orig'], para_vs['ks2_none']
    m0_inf, m0_vs, m0_va, m0_sa, m0_as, m0_orig, m0_none = para_vs['m0_inf'], para_vs['m0_vs'], para_vs['m0_va'], para_vs['m0_sa'], para_vs['m0_as'], para_vs['m0_orig'], para_vs['m0_none']
    yd_inf, yd_vs, yd_va, yd_sa, yd_as, yd_orig, yd_none = para_vs['yd_inf'], para_vs['yd_vs'], para_vs['yd_va'], para_vs['yd_sa'], para_vs['yd_as'], para_vs['yd_orig'], para_vs['yd_none']

    print ('Total scan rounds:')
    print (len(ls_novh)+len(ls_nomh)+len(ls_vh)+len(ls_inf)+len(ls_vs)+len(ls_va)+len(ls_sa)+len(ls_as)+len(ls_orig)+len(ls_none))
    # Draw tree-level boundaries
    '''
    ls_min, ls_max = None, None
    ls_lim = ls_range_IIa(lh, vh2, vs2, ks2)
    if ls_lim is not None:
        ls_min = ls_lim[0] * np.ones(100)
        ls_max = ls_lim[1] * np.ones(100)
    '''
    vspr = np.linspace(xmin, xmax, 100)
        
    # Plot args
    size = 10.
    evenly_spaced_interval = range(8)
    colors = [plt.cm.Dark2(x) for x in evenly_spaced_interval]
 
    ma_vh = []   
    for i in range(len(lsh_vh)):
	ma = (0.5*lsh_vh[i]*246**2-ls_vh[i]*vs2_vh[i]-2*ks2_vh[i])**.5
	ma_vh.append(ma)

    figma = plt.figure(figsize=(8,6))
    axma = figma.add_subplot(111)
    axma.scatter(ma_vh, lsh_vh, s=size, c='red', label='Meet BC')
    axma.legend(fontsize=10, loc='upper right', bbox_to_anchor=(1, 1))
    axma.set_xlabel(r'$m_a\ [GeV]$', size=15)
    axma.set_ylabel(r'$\lambda_{SH}$', size=15)
    plt.savefig('%s_ma_lsh.pdf' % FILE)

    #fig = plt.figure(figsize=(14,6))
    #spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2,1])
    #ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
    #ax.scatter(np.sqrt(vs2_novh+2*np.divide(ks2_novh,ls_novh)), lsh_novh, s=size, c='grey', alpha=0.5, label='Higgs vev is not a minimum')
    #ax.scatter(np.sqrt(vs2_nomh+2*np.divide(ks2_nomh,ls_nomh)), lsh_nomh, s=size, c='black', alpha=0.8, label=r'$m_h|_{(v_H,0,0)}\neq 125GeV$')
    ax.scatter(np.sqrt(vs2_inf+2*np.divide(ks2_inf,ls_inf)), lsh_inf, s=size, c=colors[0], alpha=0.5, label='Potential not bounded from below')
    ax.scatter(np.sqrt(vs2_va+2*np.divide(ks2_va,ls_va)), lsh_va, s=size, c='cyan', alpha=0.8, label=r'$V(0,0,a)<V(h,0,0), |a|>0$')
    ax.scatter(np.sqrt(vs2_vs+2*np.divide(ks2_vs,ls_vs)), lsh_vs, s=size, c=colors[2], alpha=0.8, label=r'$V(0,s,0)<V(h,0,0), |s|>0$')
    ax.scatter(np.sqrt(vs2_sa+2*np.divide(ks2_sa,ls_sa)), lsh_sa, s=size, c=colors[3], alpha=0.8, label=r'$V(0,s,a)<V(h,0,0), |s|>|a|>0$')
    ax.scatter(np.sqrt(vs2_as+2*np.divide(ks2_as,ls_as)), lsh_as, s=size, c=colors[5], alpha=0.5, label=r'$V(0,s,a)<V(h,0,0), |a|>|s|>0$')
    ax.scatter(np.sqrt(vs2_orig+2*np.divide(ks2_orig,ls_orig)), lsh_orig, s=size, c='pink', alpha=0.8, label=r'$V(0,0,0)<V(h,0,0)$')
    ax.scatter(np.sqrt(vs2_none+2*np.divide(ks2_none,ls_none)), lsh_none, s=size, c='gray', alpha=0.8, label=r'Other global minimum')
    ax.scatter(np.sqrt(vs2_vh+2*np.divide(ks2_vh,ls_vh)), lsh_vh, s=size, c='red', label='Satisfy T=0 boundary conditions')
    ax.plot(vspr, [lsh_min_bfb(vsp, m0, yd, ma) for vsp in vspr], c='blue', linewidth = 2.0, label=r'$\lambda_{SH,min}$ from BFB')
    ax.plot(vspr, [lsh_max_IIa(vsp, ma) for vsp in vspr], c='orange', linewidth = 2.0, label=r'$\lambda_{SH,max}$ from 4(a)')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(5e-3, ymax)
    ax.set_title(r'$m_a=$%sGeV' % int(ma))
    ax.set_xlabel(r'$v^{\prime}_S$ [GeV]', size=16)
    ax.set_ylabel(r'$\lambda_{SH}$', size=16)
    #plt.savefig('%s_vsp_lsh_2.pdf' % (FILE))
    
    '''
    fig2 = plt.figure(figsize=(12,6))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2,1])
    ax2 = fig2.add_subplot(spec2[0,0])
    #ax2.scatter(m0_novh, yd_novh, s=size, c='grey', alpha=0.5, label='Wrong Higgs vev')
    #ax2.scatter(m0_nomh, yd_nomh, s=size, c='black', alpha=0.5, label='Wrong Higgs mass')
    ax2.scatter(m0_inf, yd_inf, s=size, c=colors[0], alpha=0.5, label='Potential not bounded from below')
    ax2.scatter(m0_va, yd_va, s=size, c='cyan', alpha=0.5, label=r'$V(0,0,a)<V(h,0,0), |a|>0$')
    ax2.scatter(m0_vs, yd_vs, s=size, c=colors[2], alpha=0.5, label=r'$V(0,s,0)<V(h,0,0), |s|>0$')
    ax2.scatter(m0_sa, yd_sa, s=size, c=colors[3], alpha=0.7, label=r'$V(0,s,a)<V(h,0,0), |s|>|a|>0$')
    ax2.scatter(m0_as, yd_as, s=size, c=colors[5], alpha=0.7, label=r'$V(0,s,a)<V(h,0,0), |a|>|s|>0$')
    ax2.scatter(m0_orig, yd_orig, s=size, c=colors[7], alpha=0.5, label=r'$V(0,0,0)<V(h,0,0)$')
    #ax2.scatter(m0_none, yd_none, s=size, c=colors[6], alpha=0.5, label=r'Wrong global minimum')
    ax2.scatter(m0_vh, yd_vh, s=size, c='red', label='Satisfy T=0 boundary conditions')
    ax2.legend(fontsize=10, loc='upper right', bbox_to_anchor=(1.55, 1))
    ax2.set_xlabel(r'$m_0\ [GeV]$', size=15)
    ax2.set_ylabel(r'$\lambda$', size=15)
    plt.savefig('%s_m0_yd.pdf' % (FILE))
    '''

#===================================================
#Main function
#===================================================
fig = plt.figure(figsize=(14,8))
spec = gridspec.GridSpec(ncols=2, nrows=2, height_ratios=[3,1])
ax1, ax2 = fig.add_subplot(spec[0,0]), fig.add_subplot(spec[0,1])
ax_leg = fig.add_subplot(spec[1, 0:])
ax_leg.axis('off')
plot_min(FILE_1, 15, ax1, 600., 2000., 0.03)
plot_min(FILE_2, 150, ax2, 10., 700., 8.)
handles, labels = ax2.get_legend_handles_labels()
#fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.2, -0.27), ncol=4)
#fig.savefig("%s/ma15-150_vsp_lsh.pdf" % OUTDIR)
#fig_leg = plt.figure(figsize=(14,2))
#ax = fig_leg.add_subplot(111)
#ax.axis('off')
ax_leg.legend(handles, labels, bbox_to_anchor=(1.05, 1.0), ncol=3, fontsize=14)
fig.savefig("%s/ma15-150_vsp_lsh.pdf" % OUTDIR)
