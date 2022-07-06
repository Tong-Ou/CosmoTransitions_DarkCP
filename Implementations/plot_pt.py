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

FILE = sys.argv[1]
#FILE2 = sys.argv[2]
para_dict = dd.io.load('%s_plot.h5' % FILE)
para_dict_nuc = dd.io.load('%s_plot_nuc.h5' % FILE)
#para_dict_2 = dd.io.load('%s_plot.h5' % FILE2)
#para_dict_nuc_2 = dd.io.load('%s_plot_nuc.h5' % FILE2)
'''
for key in para_dict.keys():
    try:
        para_dict[key].extend(para_dict_2[key])
        para_dict_nuc[key].extend(para_dict_nuc_2[key])
    except:
        continue
'''
def lshMin(ls, ma):
    vh2 = 246.**2
    lh = 0.129
    ma2 = ma**2
    #return (ma2+np.sqrt(ma2**2 + 4*ls*lh*vh2**2))/vh2
    
    Tew = 140.
    lh = 0.129
    ch = lh*vh2/Tew**2
    ma2 = ma**2.
    Tc = 80.
    return ((2*ma2 - 1/3*ls*Tc**2) + np.sqrt((2*ma2 - 1/3*ls*Tc**2)**2 + 16*ls*vh2*ch*(Tew**2-Tc**2)))/(2*vh2)

def lshtolsMax(musp2):
    delta = 0.5
    yd = 0.4
    ma2 = 30**2
    lh = 0.129
    vh = 246.
    T = 80.
    lsh = 2*(musp2 + ma2)/(vh**2)
    ch = 0.39 + lsh/12  
    h_EW = vh*(1-ch*T**2/(2*lh*vh**2))
    ls = (lsh*(h_EW**2/2+delta*(vh**2/2-T**2/6)) - delta*(3*musp2+yd**2*T**2/4+ma2))*3/(delta*T**2)
    print (ls)
    return lsh/ls
    
def plot_fopt(para_dict, suffix, colors):

    # Load data
    #para_dict = dd.io.load('%s_plot.h5' % FILE)

    fig = plt.figure(figsize=(16,6))
    spec = gridspec.GridSpec(ncols=2, nrows=1)
    ax1 = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])

    # Plot args
    size = 10.
    size2 = 5.
    al1 = 0.5
    al2 = 0.5
    al3 = 1.
    #evenly_spaced_interval = np.linspace(0, 1, 85)
    #colors = [plt.cm.RdYlBu(x) for x in evenly_spaced_interval]
    cc = {}

    # Fixed parameters
    vh = 246.
    vh2 = vh**2.
    lh = 0.129
    vre2 = 172.9**2
    ptlist = ['a', 'bI', 'bII', 'c']	

    tags = ['nopt']
    for i in ptlist:
	tags.append(i)

    for i in ptlist:
	for j in ptlist:
	    tags.append(i+j)

    for i in ptlist:
	for j in ptlist:
	     for k in ptlist:
		tags.append(i+j+k)

    for i in ptlist:
        for j in ptlist:
             for k in ptlist:
		 for l in ptlist:
                     tags.append(i+j+k+l)

    
    for tag in tags:
	try:
	    paras = para_dict[tag]
	except:
	    continue
        if len(paras) > 0:
	    ma2, mh2, ms2, mh20 = [], [], [], []
	    pT = []
	    delta_as, vspList, musp2List = [], [], []
	    delta = []
	    
	    for para in paras:
		print('vs2:%s ls:%s lsh:%s ks2:%s yd:%s thetaY:%s m0:%s' % (para[1],para[2],para[3],para[4],para[6],para[7],para[5]))
		m = bmt.model(246.**2, para[1], 0.129, para[2], para[3], para[4], para[6], para[7], para[5], vre2)
		
		Aa = (2**0.5)*m.m0*m.yd/2
		cS = m.cs()
		cH = m.ch()
		if m.vs2 + 2*m.ks2/m.ls >= 0:
		    vsp = (m.vs2 + 2*m.ks2/m.ls)**0.5
		else:
		    vsp = np.nan
		a = (m.vh2*cH**2)/(2*125**2) + Aa**2/(4*m.lsh*m.vh2-8*m.ls*vsp**2) - (Aa-2*cS*vsp)**2/(16*m.ls*vsp**2)
		b = 0.5*(cH*m.vh2 - cS*vsp**2 + Aa*vsp)
		c = 0.25*(-m.lh*m.vh2**2 + m.ls*vsp**4)
		delta.append(c)
		delta_as.append(Aa - 2*cS*vsp)
		musp2List.append(m.ls*m.vs2 + 2*m.ks2)
		vspList.append(vsp)
		p = cS*100**2/m.ls - vsp**2
		q = Aa*100**2/(2*m.ls)
		pT.append(4*p**3 + 27*q**2)
		ma2.append(0.5*m.lsh*246**2 - m.ls*vsp**2)
		ma = (0.5*m.lsh*246**2 - m.ls*vsp**2)**.5
		print ("ma: %s" % ma)

		'''	        	
		hess = m.d2V([vh, 0., 0.], T=0.) 
		ma2.append(hess[2][2])
		ms2.append(hess[1][1])
		#print ('ma2: %s ms2: %s' % (hess[2][2], hess[1][1]))
		mh20.append(hess[0][0])
		
		try:
		    #Tc = np.sqrt((m.lh*m.vh2 - 0.5*m.lsh*vsp**2)/cH)
		    hess1 = m.d2V([0., 0., 0.], T=100)
	            pT.append((hess1[2][2])/m.ls)
		except:
		    pT.append(np.nan)
	        '''
            paras = np.array(paras)
            vs2, ls, lsh, ks2, m0, yd = paras[:, 1], paras[:, 2], paras[:, 3], paras[:, 4], paras[:, 5], paras[:, 6]
	    #mus2 = np.multiply(vs2, ls)
	    maList = np.sqrt(ma2)
	    lshtols = np.divide(lsh, ls)
	    #ks2tols = np.divide(ks2, ls)
	    try:
                p = ax1.scatter(vspList, lsh, s=size, c=colors[tag], alpha=al2, label=tag)
		#p = ax1.scatter(ms2_ma2, mus2, s=size, c=colors[tag], alpha=al2, label=tag)
		ax2.scatter(maList, lsh, s=size, c=colors[tag], alpha=al2, label=tag)
		#ax2_new.scatter(ma2, mh2, s=size2, c=color[tag], alpha=al2, label=tag)
	    except:
		if tag == 'nopt':
		    #p = ax1.scatter(ms2_ma2, mus2, s=size, c='yellow', alpha=al2, label=tag)
		    p = ax1.scatter(vspList, lsh, s=size, c='yellow', alpha=al1, label=tag)
		    ax2.scatter(maList, lsh, s=size, c='yellow', alpha=al1, label=tag)
		    #ax2_new.scatter(ma2, mh2, s=size2, c='yellow', alpha=al1, label=tag)
		elif tag == 'a':
		    p = ax1.scatter(vspList, lsh, s=size, c='red', label=tag)
		    ax2.scatter(maList, lsh, s=size, c='red', label=tag)
		else:
		    p = ax1.scatter(vspList, lsh, s=size, alpha=al2, label=tag)
		    #p = ax1.scatter(ms2_ma2, mus2, s=size, alpha=al2, label=tag)
	 	    ax2.scatter(maList, lsh, s=size, alpha=al2, label=tag)
		    #ax2_new.scatter(ma2, mh2, s=size2, alpha=al2, label=tag)
	    cc.update({tag:p.get_facecolor()[0]})

    #lsmin, lsmax = 0, 0.005
    mar = np.linspace(1, 15, 50)
    vh2 = 246**2
    #ax1.plot(mar, [2*ma**2/vh2 for ma in mar], c='black', linewidth = 2.5, label=r'$2m_a^2/v_H^2$')
    #ax1.plot(lsr, [lshMin(ls, 15.) for ls in lsr], c='grey', linewidth = 2.5, label=r'$\lambda_{SH}^{min}\ m_a=15GeV$')
    #ax1.set_ylim(0., 0.015)
    #ax1.set_xlim(-1e-4, 3e-4)
    #ax1.set_yscale('log')
    ax1.set_xlabel(r'$\sqrt{v_S^2+2\kappa_S^2/\lambda_S}$', size=15)
    ax1.set_ylabel(r'$\lambda_{SH}$', size=15)
    ax1.legend(fontsize=10, loc='upper right', bbox_to_anchor=(0.9, 1))
    #ax1.set_ylabel(r'$\mu_S^2=\lambda_S v_S^ v 2$', size=15)

    ax2.legend(fontsize=10, loc='upper right', bbox_to_anchor=(0.9, 1))
    #ax2.set_xlabel(r'$\partial^2 V_0/\partial a^2|_{[v_H, 0, 0]}$', size=15)
    #ax2.set_xlabel(r'$\sqrt{v_S^2+2\kappa_S^2/\lambda_S}$', size=15)
    ax2.set_xlabel(r'$m_a\ [GeV]$', size=15) 
    ax2.set_ylabel(r'$\lambda_{SH}$', size=15)
    #ax2.set_ylabel(r'$\partial^2 V_0/\partial h^2|_{[0, 0, v_S^\prime]}$', size=15)
    #ax2.set_xlim(0, 750)
    #ax2.set_ylim(-1e4, 12e4)
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    plt.savefig('%s_vsp_ma_lsh%s.pdf' % (FILE, suffix))
    '''
    fig2 = plt.figure(figsize=(8,6))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1)
    ax2 = fig2.add_subplot(spec2[0,0])
    
    # Aa-to-lsh/ls
    for tag in tags:
	try:
            paras = para_dict[tag]
	except:
	    continue
        if len(paras) > 0:
            paras = np.array(paras)
            ls, lsh, m0, yd = paras[:, 2], paras[:, 3], paras[:, 5], paras[:, 6]
	    Aa = []
	    lshtols = []
	    for n in range(len(paras)):
		Aa.append(m0[n]*yd[n])
		lshtols.append(lsh[n]/ls[n])
	    try:
		ax2.scatter(lshtols, Aa, s=size, c=colors[tag], alpha=al1, label=tag)
	    except:
		if tag == 'nopt':
		    ax2.scatter(lshtols, Aa, s=size, alpha=al1, c='yellow', label=tag)
		else:
                    ax2.scatter(lshtols, Aa, s=size, alpha=al1, label=tag)


    ax2.legend(fontsize=15, loc='upper right', bbox_to_anchor=(0.9, 1))
    ax2.set_xlabel(r'$\lambda_{SH}/\lambda_S$', size=15)
    ax2.set_ylabel(r'$A_a$', size=15)
    ax2.set_xlim(100, 1000)
    plt.savefig('%s_lshtols_Aa%s.pdf' % (FILE, suffix))
    '''
    return cc
    

#===================================================
#Main function
#===================================================
cc = plot_fopt(para_dict, '', {})
cc2 = plot_fopt(para_dict_nuc, '_nuc', {})

