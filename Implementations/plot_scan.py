'''
author: Tong Ou
This script is to plot the viable parameter space based on results of ran_scan.py.
'''

import numpy as np
import matplotlib.pyplot as plt
import sys 
import deepdish as dd
import argparse

ncpu = 4

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

