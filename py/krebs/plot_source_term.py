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
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os,sys
import numpy as np
import extensions

import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


if __name__ == '__main__':
  f = h5py.File('sources.h5', 'r')
  
  ds_p = f['p']
  ds_o2 = f['o2']
  ds_S = f['S']
  # axes are S, p, o2
  ds_q = f['q']
  params = f['parameters'].attrs
  kdeath = -params['death_coefficient']
  kprol = params['proliferation_coefficient']  
  p0_tum = params['homeostatic_pressure_tumor']
  p0_norm = params['homeostatic_pressure']
  pw = params['pressure_sensitivity']
  minp = min(p0_tum, p0_norm) - kprol/pw
  maxp = max(p0_tum, p0_norm) - kdeath/pw
  kpm = 0.5*(kprol - kdeath)*0.2
  ppm = 0.5*(maxp-minp)*0.2
  o2range = range(5, len(ds_o2), 5)
  
  pdfpages = PdfPages("sources.pdf")
  
  if True:
    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.set_title('proliferation rate coefficient $q(p, O_2 = 1)$')
    # full o2, no tumor
    plot.plot(ds_p, ds_q[0,:,-1], 'g', label = '$S = 0$ (normal)', linewidth=2.)
    for i in o2range:
      plot.plot(ds_p, ds_q[-1,:,i], 'g--')
    # intermediate
    plot.plot(ds_p, ds_q[len(ds_S)/2,:,-1], 'b', label = '$S = 1/2$', linewidth=2.)
    # full o2, tumor
    plot.plot(ds_p, ds_q[-1,:,-1], 'r', label = '$S = 1$ (tumor)', linewidth=2.)
    for i in o2range:
      plot.plot(ds_p, ds_q[-1,:,i], 'r--')
    # min, max
    plot.plot(ds_p, kprol*np.ones_like(ds_p), 'k')
    plot.text(0.5*(minp+maxp), kprol, r'$k^{(prol)}$')
    plot.plot(ds_p, kdeath*np.ones_like(ds_p), 'k')
    plot.text(0.5*(minp+maxp), kdeath, r'$k^{(death)}$')
    plot.plot(ds_p, np.zeros_like(ds_p), 'k', lw = 1.)
    # anotate
    plot.plot([p0_norm, p0_tum], [0., 0.], 'ko')
    plot.annotate('$p_0$ (normal)', xy=(p0_norm,0.), xycoords='data', xytext=(p0_norm - ppm, 5 * kpm), textcoords='data', arrowprops=dict(arrowstyle='->'))
    plot.annotate('$p_0$ (tumor)', xy=(p0_tum,0.), xycoords='data', xytext=(p0_tum, 3 * kpm), textcoords='data', arrowprops=dict(arrowstyle='->'))
    for i in o2range:
      o2 = ds_o2[i]
      f = 1.05
      plot.text(f*minp+(1.-f)*maxp, np.amax(ds_q[0,:,i]), r'$%s%0.1f$' % ('O_2=' if i==o2range[-1] else '', o2))
    # make it beautifull
    plot.set_ylabel('$q$')
    plot.set_xlabel('$p$')
    plot.set_xlim(minp - ppm, maxp + ppm)
    plot.set_ylim(kdeath - kpm, kprol + kpm)
    plot.grid(True)
    plot.legend()  
    pdfpages.savefig(fig)
  
  if True:
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    #plot = fig.add_subplot(111, projection='3d')
    plot = Axes3D(fig)
    plot.set_title('proliferation rate coefficient $q(p,O_2)$')
    xx = (50, -40)
    X = np.asarray(ds_p)[xx[0]:xx[1]]
    Y = np.asarray(ds_o2)
    X, Y = np.meshgrid(X, Y)
    Z = np.asarray(ds_q[0,:,:]).transpose()[:,xx[0]:xx[1]]
    p = plot.plot_wireframe(X, Y, Z, cstride=10, rstride=5, color='g', label='$S = 0$ (normal)')
    Z = np.asarray(ds_q[len(ds_S)/2,:,:]).transpose()[:,xx[0]:xx[1]]
    p = plot.plot_wireframe(X, Y, Z, cstride=10, rstride=5, color='b', label='$S = 1/2$')
    Z = np.asarray(ds_q[-1,:,:]).transpose()[:,xx[0]:xx[1]]
    p = plot.plot_wireframe(X, Y, Z, cstride=10, rstride=5, color='r', label='$S = 1$ (tumor)')
    # beautyfy
    plot.set_xlabel('$p$')
    #plot.set_xlim(minp - ppm, maxp + ppm)
    plot.set_ylabel('$O_2$')
    #plot.set_ylim(-0.05, 1.05)
    #plot.set_zlabel('$q$')
    #plot.set_zlim(kdeath - kpm, kprol + kpm)  
    plot.legend()
    pdfpages.savefig(fig)
  
  if False:
    #bad image plot
    fig = plt.figure()
    plot = fig.add_subplot(111)
    Z = np.asarray(ds_q[0,:,:]).transpose()
    plot.imshow(Z, 
                extent=(ds_p[0],ds_p[-1],ds_o2[0],ds_o2[-1]),
                origin='lower',
                interpolation='bicubic',
                cmap = matplotlib.cm.spring,
                aspect = 0.75*(maxp-minp)/(ds_o2[-1]-ds_o2[0])
                )
    X = np.asarray(ds_p)
    Y = np.asarray(ds_o2)
    X, Y = np.meshgrid(X, Y)
    levels = np.linspace(kdeath+kpm*0.1, kprol-kpm*0.1, 10)
    #[kdeath+kpm*0.1, kdeath*0.5, 0., kprol*0.5, kprol-kpm*0.1]
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    c = plot.contour(X, Y, Z
                     ,levels
                     ,colors = 'k')
    #c.collections[2].set_linewidth(2.)
    plot.clabel(c
                ,fontsize='6pt')
    plot.set_xlabel('$p$')
    plot.set_xlim(minp - ppm, maxp + ppm)
    plot.set_ylabel('$O_2$')
    plot.set_ylim(-0.01, 1.01)
    
  #plt.show()
  pdfpages.close()
  
