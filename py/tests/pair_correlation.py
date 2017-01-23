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
from os.path import basename, dirname, join, splitext
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import numpy as np
import scipy
import scipy.spatial
import scipy.fftpack
import scipy.io
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib
from third_party.utilities import *
from third_party.paircorrelation import pairCorrelationFunction_3D

import krebsutils as ku
import qsub
import h5py
import myutils
import krebs.analyzeGeneral
import krebs.plotIff

import pydicom
def std_example():
  # Particle setup
  domain_size = 20.0
  num_particles = 1000
  
  # Calculation setup
  dr = 0.1
  
  ### Random arrangement of particles ###
  particle_radius = 0.1
  rMax = domain_size / 4
  x = np.random.uniform(low=0, high=domain_size, size=num_particles)
  y = np.random.uniform(low=0, high=domain_size, size=num_particles)
  z = np.random.uniform(low=0, high=domain_size, size=num_particles)
  
  # Compute pair correlation
  g_r, r, reference_indices = pairCorrelationFunction_3D(x, y, z, domain_size, rMax, dr)
  
  # Visualize
  plt.figure()
  plt.plot(r, g_r, color='black')
  plt.xlabel('r')
  plt.ylabel('g(r)')
  plt.xlim( (0, rMax) )
  plt.ylim( (0, 1.05 * g_r.max()) )
  plt.show()
def from_vessel_file(filenames,grp_pattern):
  dirs = set()
  dataman = myutils.DataManager(20, [krebs.plotIff.DataTissue(), krebs.plotIff.DataGlobalIff(), krebs.plotIff.DataRadialIff(), krebs.analyzeGeneral.DataDistanceFromCenter(), krebs.analyzeGeneral.DataBasicVessel()])

  #run with grp_pattern: iff/vessels
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      d = myutils.walkh5(f, grp_pattern)
      assert len(d), 'you fucked up, pattern "%s" not found in "%s"!' % (grp_pattern, fn)
      dirs =set.union(dirs, d)
      for group_path in dirs:
        vesselgroup = f[group_path]
        ldvessels = ku.read_lattice_data_from_hdf(vesselgroup['lattice'])
        fieldld = ku.SetupFieldLattice(ldvessels.worldBox, 3, 10, 0.)        
        #fieldldFine = ku.SetupFieldLattice(fieldld.worldBox, 3, 50 / 10, 0.)        
        phi_vessels = krebs.analyzeGeneral.CalcPhiVessels(dataman, f['iff/vessels'], fieldld, scaling = 1., samples_per_cell = 5)        
        #a_abs = np.fabs(phi_vessels)
        #a_max = np.amax(a_abs)
        #n_phi = phi_vessels/a_max
        #visual =  (phi_vessels-phi_vessels.min())
        #inner1 = matplotlib.cm.gist_earth(n_phi)
        #inner2 = np.uint8(inner1*255)
        #result = Image.fromarray(inner2[:,:,:,0:2],mode='RGB')
        #result.save('out.tiff')
        print('bla')
        import nibabel as nib
        new_image=nib.Nifti1Image(phi_vessels,affine=np.eye(4))
        new_image.to_filename('myni')
        #for i in np.arange(0,phi_vessels.shape[2]):
        #  plt.imsave('img/%03i.png' % i,phi_vessels[:,:,i])
        #scipy.io.savemat('test.mat', phi_vessels)
        #np.save('mydata.tiff',phi_vessels)
        #pydicom.
        #pydicom.write_file('mydata.dcm', phi_vessels)      
        plt.imshow(phi_vessels[:,:,20])
        plt.show()

  
if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute po2 distributions and analyze them. In normal mode it takes arguments: parameter set name, filenames, h5 group pattern. In analysis mode -a, it takes arguments: filenames, h5 group pattern.')  
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')   
  parser.add_argument('grp_pattern',help='Where to find the vessel group in the file')  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
  

  from_vessel_file(filenames, goodArguments.grp_pattern)