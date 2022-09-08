
import os
import sys
import deepdish as dd
import baseMo_s_cpv as bmt

ntasks = int(sys.argv[1])
FILE = sys.argv[2]
#OUTDIR = FILE.split('pt')[0]
#phases = dd.io.load('%s/phases.h5' % OUTDIR)
cpv_paras = open('%s_cpvparas.dat' % FILE, 'w')
cpv_paras_Tc = open('%s_cpvparasTc.dat' % FILE, 'w')
nuc_list = open('%s_nucList.dat' % FILE, 'w')

# a: EWSFOPT, b: EWFOPT, c: EWSOPT, d: Not EWPT
tag0 = ''
para_dict = {'nopt':[]}
para_dict_nuc = {'nopt':[]}
para_lows, para_highs, para_negs = [], [], []
ptlist = ['a', 'bI', 'bII', 'c']
for i in ptlist:
  tag1 = tag0 + i
  para_dict.update({tag1:[]})
  para_dict_nuc.update({tag1:[]})
  for j in ptlist:
    tag2 = tag1 + j
    para_dict.update({tag2:[]})
    para_dict_nuc.update({tag2:[]})
    for k in ptlist:
      tag3 = tag2 + k
      para_dict.update({tag3:[]})
      para_dict_nuc.update({tag3:[]})

def categorizePt(pts):
  tag = '' 
  for pt in pts:
    hkey, lkey, hvev, lvev, Tc, pt_type = pt[0], pt[1], pt[2], pt[3], pt[4], pt[5]
    if pt_type == 2:
      tag += 'c'
    else:
      if hvev/(Tc+1e-4) < .1 and lvev/(Tc+1e-4) > 1.:
        tag += 'a'
      elif abs(hvev-lvev) < .1:
	tag += 'c'
      elif hvev/(Tc+1e-4) > .1 and lvev/(Tc+1e-4) > 1.:
        tag += 'bI'
      elif hvev/(Tc+1e-4) < .1 and lvev/(Tc+1e-4) < 1.:
	tag += 'bII'
    '''
    if lvev < 1.:
      tag += 'd'
    else:
      if pt_type == 2:
	tag += 'c'
      else:
	if hvev/Tc < .5 and lvev/Tc > 1.:
	  tag += 'a'
	else:
	  tag += 'b'
  '''
  return tag 

vs2, ls, lsh, ks2, m0, yd, thetaY = None, None, None, None, None, None, None
highPhase, zeroPhase = 0, 0
bm = None
for n in range(ntasks):
  try:
    log = open('%s_%s.log' % (FILE, n), 'r')
  except:
    continue
  lines = log.readlines()
  i = 0
  while i < len(lines):
    if 'Index' in lines[i]:
      print (lines[i])
      paras = lines[i].split()
      index = int(paras[1])
      vh2 = float(paras[2].split('vh2:')[1])
      vs2 = float(paras[3].split('vs2:')[1])
      lh = float(paras[4].split('lh:')[1])
      ls = float(paras[5].split('ls:')[1])
      lsh = float(paras[6].split('lsh:')[1])
      ks2 = float(paras[7].split('ks2:')[1])
      yd = float(paras[8].split('yd:')[1])
      thetaY = float(paras[9].split('thetaY:')[1])
      m0 = float(paras[10].split('m0:')[1])
      v2re = float(paras[11].split('v2re:')[1])
      bm = bmt.model(vh2, vs2, lh, ls, lsh, ks2, yd, thetaY, m0, v2re)
      #bm.phases = phases[index]
      pts = []
      pts_nuc = []
      high_keys, low_keys = [], []
      high_keys_nuc, low_keys_nuc = [], []
      action, delta_rho = None, None

    if 'High-T phase' in lines[i]:
      try:
        highPhase = int(lines[i].split()[-1])
	print ('High Phase: %s' % highPhase)
      except:
        pass

    if 'Zero-T phase' in lines[i]:
      try:
        zeroPhase = int(lines[i].split()[-1])
	print ('Zero Phase: %s' % zeroPhase)
      except:
        pass
      
    if 'transition at Tc' in lines[i]:
      Tc = float(lines[i].split('Tc = ')[1])
      high_line = lines[i+2].split()
      try:
        high_key = int(high_line[2].split(';')[0])
      except:
	i += 7
	continue
      try:
        high_hvev = float(high_line[5].split('[')[1])
	high_svev = float(high_line[6])
	high_avev = float(high_line[7].split(']')[0])
      except:
	high_hvev = float(high_line[6])
	high_svev = float(high_line[7])
	high_avev = float(high_line[8].split(']')[0])
      low_line = lines[i+4].split()
      try:
        low_key = int(low_line[2].split(';')[0])
      except:
	i+=7
	continue
      try:
        low_hvev = float(low_line[5].split('[')[1])
	low_svev = float(low_line[6])
	low_avev = float(low_line[7].split(']')[0])
      except:
	low_hvev = float(low_line[6])
	low_svev = float(low_line[7])
	low_avev = float(low_line[8].split(']')[0])
      delta_rho = float(lines[i+5].split()[3])

      #pts.append([high_key, low_key, high_hvev, low_hvev, Tc])
      # If this transition can happen (do not consider nucleation)
      if high_key == highPhase or high_key in low_keys:
        if not(high_key in high_keys):
          high_keys.append(high_key)
	  low_keys.append(low_key)
	  if 'First' in lines[i]:
	    pt_type = 1
	  else:
	    pt_type = 2
	  pts.append([high_key, low_key, high_hvev, low_hvev, Tc, pt_type, high_svev, high_avev, low_svev, low_avev, delta_rho])
          print ('Found transition from %s to %s at Tc %s' % (high_key, low_key, Tc)) 
	  print ('\n')
      i += 7
      continue

    if 'The T-dependent phase potential reads' in lines[i]:
      paras = [index, vs2, ls, lsh, ks2, m0, yd, thetaY]
      if len(low_keys) == 0 or not (low_keys[-1] == zeroPhase):
	print ('nopt\n')
	para_dict['nopt'].append(paras)
      else:
	tag = categorizePt(pts)
	try:
	  para_dict[tag].append(paras)
	except:
	  para_dict.update({tag:[paras]})
	if tag[-1] == 'a':
	  paras.append(pts[-1][6])
	  paras.append(pts[-1][7])
	  paras.append(pts[-1][4])
	  paras.append(pts[-1][10])
	  for p in paras:
	    cpv_paras_Tc.write('%s ' % p)
	  cpv_paras_Tc.write('\n')
	 
	print (tag+'\n')

    # Consider nucleation
    if 'transition at Tnuc' in lines[i]:
      Tnuc = float(lines[i].split()[-1])
      high_line = lines[i+2].split()
      try:
        high_key = int(high_line[2].split(';')[0])
      except:
        i += 10
        continue
      try:
        high_hvev = float(high_line[5].split('[')[1])
	high_svev = float(high_line[6])
	high_avev = float(high_line[7].split(']')[0])
      except:
        high_hvev = float(high_line[6])
	high_svev = float(high_line[7])
	high_avev = float(high_line[8].split(']')[0])
      low_line = lines[i+4].split()
      try:
        low_key = int(low_line[2].split(';')[0])
      except:
        i += 10
        continue
      try:
        low_hvev = float(low_line[5].split('[')[1])
	low_svev = float(low_line[6])
	low_avev = float(low_line[7].split(']')[0])
      except:
        low_hvev = float(low_line[6])
	low_svev = float(low_line[7])
	low_avev = float(low_line[8].split(']')[0])
      delta_rho = float(lines[i+6].split()[3])
      action = float(lines[i+8].split()[-1])
      #dsdt = 1e100
      dsdt = float(lines[i+9].split()[-1])

      if high_key == highPhase or high_key in low_keys_nuc:
        if all((not(high_key in high_keys_nuc), action <= 150., action >= 130.)):
          high_keys_nuc.append(high_key)
	  low_keys_nuc.append(low_key)
	  if 'First' in lines[i]:
	    pt_type = 1
	  else:
	    pt_type = 2

	  xlow = [low_hvev, low_svev, low_avev]
	  xhigh = [high_hvev, high_svev, high_avev]
	  #dsdt = bm.dSdT(xlow, xhigh, Tnuc)
          pts_nuc.append([high_key, low_key, high_hvev, low_hvev, Tnuc, pt_type, high_svev, high_avev, low_svev, low_avev, delta_rho])
	  #pts_nuc.append([high_key, low_key, high_hvev, low_hvev, Tnuc, pt_type, high_svev, high_avev, low_svev, low_avev, delta_rho])
	  nuc_list.write('%s ' % index)
	  print ('Found nucleation between %s and %s at Tnuc=%s' % (high_key, low_key,Tnuc))
	  print ('Vev of s before phase transition: %s' % (high_svev))
	  print ('Vev of a before phase transition: %s' % (high_avev))
	  print ('\n')
      i += 10
      continue

    # Make sure the parameters are recorded only once
    if any(('The parameters are:' in lines[i] and i>0, i == len(lines)-1)):
      paras = [index, vs2, ls, lsh, ks2, m0, yd, thetaY]
      if len(low_keys_nuc) == 0 or not(low_keys_nuc[-1] == zeroPhase) :
	if action is None:
	  para_dict_nuc['nopt'].append(paras)
	elif all((action > 0., action < 130.)):
	  para_lows.append(paras)
	elif action > 150.:
	  para_highs.append(paras)
	elif action < 0.:
	  para_negs.append(paras)
      else: 	
	tag = categorizePt(pts_nuc)
	if tag.endswith('a'):
	  high_s, high_a, high_h, low_s, low_a, low_h, Tnuc, delta_rho = pts_nuc[-1][6], pts_nuc[-1][7], pts_nuc[-1][2], pts_nuc[-1][8], pts_nuc[-1][9], pts_nuc[-1][3], pts_nuc[-1][4], pts_nuc[-1][10]
	  #high_s, high_a, Tnuc, delta_rho = pts_nuc[-1][6], pts_nuc[-1][7], pts_nuc[-1][4], pts_nuc[-1][10]
	  paras.append(high_s)
	  paras.append(high_a)
	  paras.append(high_h)
	  paras.append(low_s)
	  paras.append(low_a)
	  paras.append(low_h)
	  paras.append(Tnuc)
	  paras.append(delta_rho)
	  paras.append(dsdt)
	  for p in paras:
	    cpv_paras.write('%s ' % p)
	  cpv_paras.write('\n')
	para_dict_nuc[tag].append(paras)
    i += 1    

para_dict_nuc.update({'lows':para_lows})
para_dict_nuc.update({'highs':para_highs})
para_dict_nuc.update({'negs':para_negs})  

filename = '%s_plot.h5' % FILE
filename_nuc = '%s_plot_nuc.h5' % FILE
dd.io.save(filename, para_dict)
dd.io.save(filename_nuc, para_dict_nuc)
cpv_paras.close()
cpv_paras_Tc.close()
nuc_list.close()
