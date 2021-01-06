'''

author: Tong Ou
This script is to study the upper bound for the
CP violation factor ks2.

'''
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

FILE = sys.argv[1]
#OUT_PATH = sys.argv[2]

ncpu = 4
para_list = []
for rank in range(ncpu):
    filename = '%s_%s.npy' % (FILE, rank)
    if os.path.exists(filename):
        para = np.load(filename, allow_pickle = True)
        para_list.append(para)
    else:
        continue

paras = np.concatenate(para_list,axis = 0)
np.save('%s.npy' % FILE, paras)

ks2 = []
ls = []
for para in paras:
    ks2.append(para[5])
    ls.append(para[3])

# Plotting
plt.figure()
ax = plt.subplot(111)
ax.scatter(ks2, ls, label='pass BC')
ax.legend()
ax.set_xlabel(r'$\kappa_S^2$')
ax.set_ylabel(r'$\lambda_{S}$')
ax.grid(True)
plt.show()
