#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  
import os, sys
import h5files
from krebsutils import typelist

def createParamFileFromSlurm(label):
  writepath = '/localdisk/thierry/tmp/%s_parms.py' % label
  #mode = 'a' if os.path.exists(writepath) else 'w'
  mode = 'w'
  if os.path.exists(writepath):
    os.remove(writepath)
  #get slurm files
  slurm_files=[]
  for fn in os.listdir('.'):
    if "slurm" in fn:
      slurm_files.append(fn)
  
  with open(writepath,mode) as f:
    
    f.write('# -*- coding: utf-8 -*-\n')
    for t in typelist:
      for slurm_file in slurm_files:
        with open(slurm_file,'r') as f_read:
          lines = f_read.readlines()
          if t in lines[1]:
            for (i,line) in enumerate(lines):
              #print line
              if 'Best' in line:
                list_of_words = line.split()
                list_of_words_next_line = lines[i+1].split()
                best_k_m = list_of_words[5]
                best_k_c = list_of_words[6]
                best_k_s = list_of_words[7]
                best_l = list_of_words[8][:-1]
                #best_q = list_of_words_next_line[0][:-1]
                print best_k_m,best_k_c,best_k_s,best_l, #best_q
                #for word in line.split():
                #  print word
  #          f_opt_data = h5files.open(pso_file, 'r', search = False)
  #          xopt = f_opt_data.attrs.get('xopt')
            f.write('pso_param_%s_%s_vary1=deepcopy(adaption_master_phase)\n' % (label,t))
            f.write('pso_param_%s_%s_vary1[\'adaption\'].update(\n' % (label,t))
            f.write('  k_m = %f,\n' % float(best_k_m))
            f.write('  k_c = %f,\n' % float(best_k_c))
            f.write('  k_s = %f,\n' % float(best_k_s))
            f.write('  cond_length = %f,\n' % float(best_l))
            #f.write('  Q_refdot = %f,\n' % float(best_q))
            f.write(')\n')
            print("end createParamFile")
        
def createParamFile(files, label):
  writepath = '/localdisk/thierry/tmp/%s_parms.py' % label
  #mode = 'a' if os.path.exists(writepath) else 'w'
  mode = 'w'
  if os.path.exists(writepath):
    os.remove(writepath)
  with open(writepath,mode) as f:
    f.write('# -*- coding: utf-8 -*-\n')
    for t in typelist:
      for pso_file in files:
        if t in pso_file:
          f_opt_data = h5files.open(pso_file, 'r', search = False)
          xopt = f_opt_data.attrs.get('xopt')
          f.write('pso_param_%s_%s_vary1=deepcopy(adaption_master)\n' % (label,t))
          f.write('pso_param_%s_%s_vary1[\'adaption\'].update(\n' % (label,t))
          f.write('  k_m = %f,\n' % xopt[0])
          f.write('  k_c = %f,\n' % xopt[1])
          f.write('  k_s = %f,\n' % xopt[2])
          f.write('  cond_length = %f,\n' % xopt[3])
          f.write('  Q_refdot = %f,\n' % xopt[4])
          f.write(')\n')
  print("Created file: %s"%writepath)
  
#p3d_mini_typeG_vary1= deepcopy(appReal21b)
#  p3d_mini_typeG_vary1['adaption'].update(
#    k_c = 2.5,#2.85,
#    k_m = 0.8,
#    k_s = 1.3,
#    Q_refdot = 40.0,
#    S_0 = 20.,
#    cond_length = 2500.,
#    #if this is 0 we iterate until qdev is reached
#    max_nun_iterations = 300,
#    qdev = .0,
#    #if starting_radii is 0. we use the values given in
#    #the input file
#    starting_radii = 0.,
#    delta_t = 0.,
#)

if __name__ == '__main__':
  import glob
  import qsub
  import myutils
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead

  sys.argv = qsub.parse_args(sys.argv)

  parser = optparse.OptionParser()
#  parser.add_option("-d","--data", dest="datalist", help="which data (pressure, flow, shearforce, hematocrit) as comma separated list", default='pressure', action="store")
  parser.add_option("-s","--slurm-files", dest="from_slurm_files", help="creates File from slurm files", default=False, action="store_true")
#  parser.add_option("--filter-radius-high-pass", dest="filterradiushighpass", action="store", type="float", default = -1)
#  parser.add_option("--filter-radius-low-pass", dest="filterradiuslowpass", action="store", type="float", default = -1)      
#  parser.add_option("--no-overlay", dest="overlay", default = True, action="store_false")
#  parser.add_option("--dpi", dest="dpi", default=None, action="store")
#  parser.add_option("--format", dest="format", default=None, action="store")
#  parser.add_option("","--cam", dest="cam", help="camera mode: topdown, pie, topdown_slice", default="topdown_slice", action="store", type="string")
  options, args = parser.parse_args()

  if(options.from_slurm_files):
    label = args[0]
    createParamFileFromSlurm(label)
  else:
    #files = args[:-1]
    label = args[-1]
    postfix = ''
    files=[]
    for fn in os.listdir('.'):
      if ".h5" in fn:
        files.append(fn)
    createParamFile(files, label)
  
  
