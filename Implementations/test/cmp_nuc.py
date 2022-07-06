import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

FILE = str(sys.argv[1])
nucFile_full = open('%s/full_potential/pt_nucList.dat' % FILE ,'r')
nucList_full = nucFile_full.readlines()[0].split()
yfull = np.ones(len(nucList_full)) # for plotting

nucFile_highT = open('%s/highT/pt_nucList.dat' % FILE, 'r')
nucList_highT = nucFile_highT.readlines()[0].split()
yhighT = np.ones(len(nucList_highT))*2

fig = plt.figure(figsize=(10,8))
spec = gridspec.GridSpec(ncols=1, nrows=2)
ax1 = fig.add_subplot(spec[0,0])
ax2 = fig.add_subplot(spec[1,0])

ax1.scatter(nucList_full, yfull, s=5, label='Full potential')
#ax1.set_xlim(0, 1000)
ax1.set_yticks([])
#ax1.set_xticks([])
ax1.legend()

ax2.scatter(nucList_highT, yhighT, s=5, c='orange', label='High T')
#ax2.set_xlim(0, 1000)
ax2.set_yticks([])
#ax2.set_xticks([])
ax2.legend()

plt.savefig('%s/nuc_cmp.pdf' % FILE)

print ('Nucleations found by full potential but not by high T:')
for nuc in nucList_full:
    if not (nuc in nucList_highT):
	print (nuc)
