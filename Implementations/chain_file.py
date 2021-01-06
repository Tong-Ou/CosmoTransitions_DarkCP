'''
This script is to chain the files produced by multiple 
CPU cores.
'''
import numpy as np
import os
import sys

ncpu = int(sys.argv[1])
FILE = sys.argv[2]
para_list = []
for i in range(ncpu):
    filename = '%s_%s.npy' % (FILE, i)
    if os.path.exists(filename):
        para = np.load(filename, allow_pickle = True)
        para_list.append(para)
    else:
        continue
paras = np.concatenate(para_list,axis = 0)
np.save('%s.npy' % FILE, paras)

