'''
This script is to chain the phases h5 files.
(for older scannings)
'''
import numpy as np
import os
import sys
import deepdish as dd

ncpu = int(sys.argv[1])
OUTDIR = sys.argv[2]

phases_dict_0 = dd.io.load('%s/phases_0.h5' % OUTDIR)
for i in range(1, ncpu):
    try:
        phases_dict = dd.io.load('%s/phases_%s.h5' % (OUTDIR, i))
        phases_dict_0.update(phases_dict)
    except:
	continue

dd.io.save('%s/phases.h5' % OUTDIR, phases_dict_0)
