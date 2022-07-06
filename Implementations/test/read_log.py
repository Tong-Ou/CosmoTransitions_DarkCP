
import os
import sys
import deepdish as dd

ntasks = int(sys.argv[1])
FILE = sys.argv[2]

ls_inf, lsh_inf, vs2_inf, ks2_inf, m0_inf, yd_inf, thetaY_inf = [], [], [], [], [], [], []
ls_vs, lsh_vs, vs2_vs, ks2_vs, m0_vs, yd_vs, thetaY_vs = [], [], [], [], [], [], []
ls_va, lsh_va, vs2_va, ks2_va, m0_va, yd_va, thetaY_va = [], [], [], [], [], [], []
ls_sa, lsh_sa, vs2_sa, ks2_sa, m0_sa, yd_sa, thetaY_sa = [], [], [], [], [], [], []
ls_as, lsh_as, vs2_as, ks2_as, m0_as, yd_as, thetaY_as = [], [], [], [], [], [], []
ls_orig, lsh_orig, vs2_orig, ks2_orig, m0_orig, yd_orig, thetaY_orig = [], [], [], [], [], [], []
ls_none, lsh_none, vs2_none, ks2_none, m0_none, yd_none, thetaY_none = [], [], [], [], [], [], []

for n in range(ntasks):
  print ('Reading log %s' % n)
  try:
    log = open('%s_%s.log' % (FILE, n), 'r')
  except:
    continue
  lines = log.readlines()
  for i in range(len(lines)-1):
    if any((i > 0 and lines[i] == '\n' and lines[i+1] == '\n', i==len(lines)-1)):
      if 'Local minimum' in lines[i-2]:
	newline = lines[i-2].split()
	h = newline[2].split('[')[1]
	try:
	  h = abs(float(newline[2].split('[')[1]))
	  s = abs(float(newline[3]))
	  a = abs(float(newline[4].split(']')[0]))
	except:
	  h = abs(float(newline[3]))
	  s = abs(float(newline[4]))
	  a = abs(float(newline[5].split(']')[0]))
	print ('\n')
	print ('Local minimum deeper than Higg minimum: %s %s %s' % (h, s, a))
	paras = lines[i-7].split()
	ls = float(paras[0].split(':')[1])
	lsh = float(paras[1].split(':')[1])
	vs2 = float(paras[2].split(':')[1])
	ks2 = float(paras[3].split(':')[1])
	m0 = float(paras[4].split(':')[1])
	yd = float(paras[5].split(':')[1])
	thetaY = float(paras[6].split(':')[1])
	print ('ls:%s lsh:%s vs2:%s ks2:%s m0:%s yd:%s thetaY:%s' % (ls, lsh, vs2, ks2, m0, yd, thetaY))
	if any((h > 1e10, s > 1e10, a > 1e10)):
	  ls_inf.append(ls)
	  lsh_inf.append(lsh)
	  vs2_inf.append(vs2)
	  ks2_inf.append(ks2)
	  m0_inf.append(m0)
	  yd_inf.append(yd)
	  thetaY_inf.append(thetaY)
	elif all((s/h > 1e3, s/a > 1e3)):
	  ls_vs.append(ls)
	  lsh_vs.append(lsh)
	  vs2_vs.append(vs2)
	  ks2_vs.append(ks2)
	  m0_vs.append(m0)
	  yd_vs.append(yd)
	  thetaY_vs.append(thetaY)
	elif all((a/h > 1e3, a/s > 1e3)):
	  ls_va.append(ls)
	  lsh_va.append(lsh)
	  vs2_va.append(vs2)
	  ks2_va.append(ks2)
	  m0_va.append(m0)
	  yd_va.append(yd)
	  thetaY_va.append(thetaY)
 	elif all((h < 1e-3, s > 0.1, a > 0.1, s > a)):
	  ls_sa.append(ls)
	  lsh_sa.append(lsh)
	  vs2_sa.append(vs2)
	  ks2_sa.append(ks2)
	  m0_sa.append(m0)
	  yd_sa.append(yd)
	  thetaY_sa.append(thetaY)
	elif all((h < 1e-3, s > 0.1, a > 0.1, s < a)):
	  ls_as.append(ls)
	  lsh_as.append(lsh)
	  vs2_as.append(vs2)
	  ks2_as.append(ks2)
	  m0_as.append(m0)
	  yd_as.append(yd)
	  thetaY_as.append(thetaY)
	elif all((h < 1e-3, s < 1e-3, a < 1e-3)):
	  ls_orig.append(ls)
	  lsh_orig.append(lsh)
	  vs2_orig.append(vs2)
	  ks2_orig.append(ks2)
	  m0_orig.append(m0)
	  yd_orig.append(yd)
	  thetaY_orig.append(thetaY)
	else:
	  print ('Local minimum does not match any known categories!')
	  ls_none.append(ls)
	  lsh_none.append(lsh)
	  vs2_none.append(vs2)
	  ks2_none.append(ks2)
	  m0_none.append(m0)
	  yd_none.append(yd)
	  thetaY_none.append(thetaY)

para_dict = {'ls_inf':ls_inf, 'lsh_inf':lsh_inf, 'vs2_inf':vs2_inf, 'ks2_inf':ks2_inf, 'm0_inf':m0_inf, 'yd_inf':yd_inf, 'thetaY_inf':thetaY_inf,
	 'ls_vs':ls_vs, 'lsh_vs':lsh_vs, 'vs2_vs':vs2_vs, 'ks2_vs':ks2_vs, 'm0_vs':m0_vs, 'yd_vs':yd_vs, 'thetaY_vs':thetaY_vs,
	 'ls_va':ls_va, 'lsh_va':lsh_va, 'vs2_va':vs2_va, 'ks2_va':ks2_va, 'm0_va':m0_va, 'yd_va':yd_va, 'thetaY_va':thetaY_va,
	 'ls_sa':ls_sa, 'lsh_sa':lsh_sa, 'vs2_sa':vs2_sa, 'ks2_sa':ks2_sa, 'm0_sa':m0_sa, 'yd_sa':yd_sa, 'thetaY_sa':thetaY_sa,
	 'ls_as':ls_as, 'lsh_as':lsh_as, 'vs2_as':vs2_as, 'ks2_as':ks2_as, 'm0_as':m0_as, 'yd_as':yd_as, 'thetaY_as':thetaY_as,
	 'ls_orig':ls_orig, 'lsh_orig':lsh_orig, 'vs2_orig':vs2_orig, 'ks2_orig':ks2_orig, 'm0_orig':m0_orig, 'yd_orig':yd_orig, 'thetaY_orig':thetaY_orig,
	 'ls_none':ls_none, 'lsh_none':lsh_none, 'vs2_none':vs2_none, 'ks2_none':ks2_none, 'm0_none':m0_none, 'yd_none':yd_none, 'thetaY_none':thetaY_none}
filename = '%s_vs.h5' % FILE
dd.io.save(filename, para_dict)
