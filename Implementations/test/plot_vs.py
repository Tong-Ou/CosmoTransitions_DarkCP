'''

author: Tong Ou
This script is to study the upper bound for the parameter vs.

'''
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import sys
import os

FILE = sys.argv[1]
#OUT_PATH = sys.argv[2]

paras = np.load('%s.npy' % FILE, allow_pickle=True)

vs = []
ls = []
for para in paras:
    vs.append(para[1]**0.5)
    ls.append(para[3])

# Plotting
plt.figure()
ax = plt.subplot(111)
ax.scatter(vs, ls, s=10, c='red', label='meet T=0 BC')
ax.legend()
ax.set_xlabel(r'$v_S$')
ax.set_ylabel(r'$\lambda_{S}$')
#ax.set_ylim(0., .2)
#ax.set_xlim(200., 500.)
ax.grid(True)
plt.savefig('%s_vs_ls.pdf' % FILE)
