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

FILE = sys.argv[1]
para_dict = dd.io.load('%s_plot.h5' % FILE)
para_dict_nuc = dd.io.load('%s_plot_nuc.h5' % FILE)

def plot_fopt(para_dict, suffix, colors):

    # Load data
    #para_dict = dd.io.load('%s_plot.h5' % FILE)

    fig = plt.figure(figsize=(8,6))
    spec = gridspec.GridSpec(ncols=1, nrows=1)
    ax1 = fig.add_subplot(spec[0,0])

    # Plot args
    size = 10.
    al1 = 0.5
    al2 = 0.8
    al3 = 1.
    #evenly_spaced_interval = np.linspace(0, 1, 85)
    #colors = [plt.cm.RdYlBu(x) for x in evenly_spaced_interval]
    cc = {}

    # One-step pt
    for tag in ['a', 'b', 'c', 'd']:
	paras = para_dict[tag]
        if len(paras) > 0:
            paras = np.array(paras)
            ls, lsh = paras[:, 1], paras[:, 2]
	    try:
                p = ax1.scatter(lsh, ls, s=size, c=colors[tag], alpha=al1, label=tag)
	    except:
		p = ax1.scatter(lsh, ls, s=size, alpha=al1, label=tag)
	    cc.update({tag:p.get_facecolor()[0]})

    # Two-step pt
    for i in ['a', 'b', 'c', 'd']:
	for j in ['a', 'b', 'c', 'd']:
	    tag = i+j
	    paras = para_dict[tag]
            if len(paras) > 0:
                paras = np.array(paras)
                ls, lsh = paras[:, 1], paras[:, 2]
		try:
		    p = ax1.scatter(lsh, ls, s=size, c=colors[tag], alpha=al2, label=tag)
		except:
                    p = ax1.scatter(lsh, ls, s=size, alpha=al2, label=tag)
		cc.update({tag:p.get_facecolor()[0]})
    
    # Three-step pt
    for i in ['a', 'b', 'c', 'd']:
	for j in ['a', 'b', 'c', 'd']:
	    for k in ['a', 'b', 'c', 'd']:
		tag = i+j+k
		paras = para_dict[tag]
                if len(paras) > 0:
            	    paras = np.array(paras)
                    ls, lsh = paras[:, 1], paras[:, 2]
		    try:
			p = ax1.scatter(lsh, ls, s=size, c=colors[tag], alpha=al3, label=tag)
		    except:
            	        p = ax1.scatter(lsh, ls, s=size, alpha=al3, label=tag)
		    cc.update({tag:p.get_facecolor()[0]})

    tag = 'nopt'
    paras = para_dict[tag]
    if len(paras) > 0:
	paras = np.array(paras)
	ls, lsh = paras[:, 1], paras[:, 2]
	try:
	    p = ax1.scatter(lsh, ls, s=size, c=colors[tag], alpha=al, label=tag)
	except:
	    p = ax1.scatter(lsh, ls, s=size, c='yellow', alpha=al1, label=tag)
	cc.update({tag:p.get_facecolor()[0]})

    ax1.legend(fontsize=10, loc='upper right', bbox_to_anchor=(0.9, 1))
    #ax1.set_ylim(0., ylim)
    #ax1.set_xlim(0., xlim)
    ax1.set_xlabel(r'$\lambda_{SH}$', size=15)
    ax1.set_ylabel(r'$\lambda_S$', size=15)
    plt.savefig('%s_lsh_ls%s.pdf' % (FILE, suffix))
    
    fig2 = plt.figure(figsize=(8,6))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1)
    ax2 = fig2.add_subplot(spec2[0,0])

    # One-step pt
    for tag in ['a', 'b', 'c', 'd']:
        paras = para_dict[tag]
        if len(paras) > 0:
            paras = np.array(paras)
            m0, yd = paras[:, 3], paras[:, 4]
	    try:
		ax2.scatter(m0, yd, s=size, c=colors[tag], alpha=al1, label=tag)
	    except:
                ax2.scatter(m0, yd, s=size, alpha=al1, label=tag)

    # Two-step pt
    for i in ['a', 'b', 'c', 'd']:
        for j in ['a', 'b', 'c', 'd']:  
            tag = i+j
            paras = para_dict[tag]      
            if len(paras) > 0:
                paras = np.array(paras)
                m0, yd = paras[:, 3], paras[:, 4]
		try:
		    ax2.scatter(m0, yd, s=size, alpha=al2, label=tag, c=colors[tag])
		except:
                    ax2.scatter(m0, yd, s=size, alpha=al2, label=tag)

    # Three-step pt
    for i in ['a', 'b', 'c', 'd']:
        for j in ['a', 'b', 'c', 'd']:
            for k in ['a', 'b', 'c', 'd']:
                tag = i+j+k
                paras = para_dict[tag]
                if len(paras) > 0:
                    paras = np.array(paras)
                    m0, yd = paras[:, 3], paras[:, 4]
		    try:
			ax2.scatter(m0, yd, s=size, alpha=al3, label=tag, c=colors[tag])
		    except:
                        ax2.scatter(m0, yd, s=size, alpha=al3, label=tag)

    tag = 'nopt'
    paras = para_dict[tag]
    if len(paras) > 0:
	paras = np.array(paras)
	m0, yd = paras[:, 3], paras[:, 4]
	try:
	    ax2.scatter(m0, yd, s=size, alpha=al1, label=tag, c=colors[tag])
	except:
	    ax2.scatter(m0, yd, s=size, alpha=al1, label=tag, c='yellow')

    ax2.legend(fontsize=10, loc='upper right', bbox_to_anchor=(0.9, 1))
    ax2.set_xlabel(r'$m_0\ [GeV]$', size=15)
    ax2.set_ylabel(r'$\lambda$', size=15)
    plt.savefig('%s_m0_yd%s.pdf' % (FILE, suffix))

    return cc
    

#===================================================
#Main function
#===================================================
cc = plot_fopt(para_dict, '', {})
cc2 = plot_fopt(para_dict_nuc, '_nuc', cc)

