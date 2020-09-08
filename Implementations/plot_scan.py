'''
author: Tong Ou
This script is to plot the viable parameter space based on results of ran_scan.py.
'''

import numpy as np
import matplotlib.pyplot as plt

ncpu = 4

ls_scan = np.zeros(0)
ls_sfo = np.zeros(0)
lsh_scan = np.zeros(0)
lsh_sfo = np.zeros(0)
sh_fo = np.zeros(0)
yd_fo = np.zeros(0)

for i in range(ncpu):
    ls_scan_i = np.load('outputs/ls_scan_%s.npy' % i)
    ls_sfo_i = np.load('outputs/ls_sfo_%s.npy' % i)
    lsh_scan_i = np.load('outputs/lsh_scan_%s.npy' % i)
    lsh_sfo_i = np.load('outputs/lsh_sfo_%s.npy' % i)
    sh_fo_i = np.load('outputs/sh_fo_%s.npy' % i)
    yd_fo_i = np.load('outputs/yd_fo_%s.npy' % i)

    if i == 0:
        ls_scan = ls_scan_i
        ls_sfo = ls_sfo_i
        lsh_scan = lsh_scan_i
        lsh_sfo = lsh_sfo_i
        sh_fo = sh_fo_i
        yd_fo = yd_fo_i
    else:
        ls_scan = np.concatenate((ls_scan, ls_scan_i), axis=0)
        ls_sfo = np.concatenate((ls_sfo, ls_sfo_i), axis=0)
        lsh_scan = np.concatenate((lsh_scan, lsh_scan_i), axis=0)
        lsh_sfo = np.concatenate((lsh_sfo, lsh_sfo_i), axis=0)
        sh_fo = np.concatenate((sh_fo, sh_fo_i), axis=0)
        yd_fo = np.concatenate((yd_fo, yd_fo_i), axis=0)

plt.figure(figsize=(12,5))
ax = plt.subplot(121)
ax.scatter(ls_scan, lsh_scan, c='black', label='scan points')
ax.scatter(ls_sfo, lsh_sfo, c='red', label=r'$v_c/T_c>1$')
ax.legend()
ax.set_xlabel(r'$\lambda_S$')
ax.set_ylabel(r'$\lambda_{SH}$')
ax.grid(True)

ax = plt.subplot(122)
ax.scatter(yd_fo, sh_fo, c='red', label='First-order Transition Points')
ax.legend()
ax.set_xlabel(r'$y_d$')
ax.set_ylabel(r'$v_c/T_c$')
ax.grid(True)

plt.show()    


