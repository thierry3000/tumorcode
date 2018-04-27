# -*- coding: utf-8 -*-

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_min_max_o2(goodArguments, pp):
  mins = []
  maxs = []
  
  print("filename: %s" % str(goodArguments.vbl_simulation_output_filename))
  with h5py.File(str(goodArguments.vbl_simulation_output_filename), 'r') as f:
    goodKeys = [str(x) for x in f.keys() if 'out' in x]
    for key in goodKeys:
      po2Field = np.asarray(f[key+"/po2/po2field"])
      mins.append(np.min(po2Field))
      maxs.append(np.max(po2Field))
      
  fig, ax = plt.subplots()    
  ax.plot(mins,label='min po2')
  ax.plot(maxs,label='max po2')
  
  
  ax.set_xlabel('#out- group')
  ax.set_ylabel('po2 value/ mmHg')
  
    
  ax.set_title('min and max value of the po2Field')
  legend = ax.legend(loc='center', shadow=True)
  if interactive:
    plt.show()
  else:
    pp.savefig()
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Plot min and max of po2field over time.')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='safe.h5', help='output file name in hdf5 format')
  
  interactive = False;
  goodArguments, otherArguments = parser.parse_known_args()
  with PdfPages('min_max_po2_Field_of_%s.pdf' % str(goodArguments.vbl_simulation_output_filename)) as pp:
    #plot_runtime_from_h5(goodArguments,pp)
    plot_min_max_o2(goodArguments,pp)