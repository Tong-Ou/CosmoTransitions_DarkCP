
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

import numpy as np

import baseMo_s_cpv as bm

l1 = 0.129
vh2 = 246.**2.
vs2 = 100.**2.
l2 = 0.1
lm = 1.6
ks2 = 3000.
yd = 0.
thetaY = 0.
m0 = 0.
v2re = 600.**2.
m = bm.model(vh2, vs2, l1, l2, lm, ks2, yd, thetaY, m0, v2re)

h = np.linspace(0.,1e10,200)
V0 = []
for x in h:
  V0.append(m.Vtot([x,0.,0.], T=0.))

plt.plot(h, V0)
#plt.xscale('log')
#plt.yscale('log')
plt.savefig('outputs/full_potential_cpv/ks23000_vs100/V0_h_lsh1p6.pdf')
