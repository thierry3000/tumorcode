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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from os.path import basename, splitext
import h5py
import os,sys
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils

import plotBulkTissue2d

import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils

f2l = lambda x: myutils.f2s(x, prec=2, latex=True)



def analyzeO2_with_tumor(f):
    
    dataman = myutils.DataManager(100, [plotBulkTissue2d.DataTumorTissueSingle()])

    
    resultfile = plotBulkTissue2d.ResultFile(f, dataman)
    group = resultfile.groupnames[-1]
    ld = resultfile.obtain_data('ld')

    print group

    arterial_reference = 12.

    theta_tumor = resultfile.obtain_data('fieldvariable', 'theta_tumor', group)
    o2 = arterial_reference * resultfile.obtain_data('fieldvariable','oxy', group)
    dist_tumor = resultfile.obtain_data('fieldvariable', 'dist_tumor', group)

    def tumorAndTissueAverages(field):
      tumor = field[theta_tumor > 0.5]
      normal = field[theta_tumor < 0.5]
      tissueavg = (np.average(normal), np.std(normal))
      tumoravg = (np.average(tumor), np.std(tumor))
      return tissueavg, tumoravg

    cons_rate_normal = float(f['parameters'].attrs['o2_cons_coeff_normal'])
    cons_rate_tumor  = float(f['parameters'].attrs['o2_cons_coeff_tumor'])
    consumption = (theta_tumor * cons_rate_tumor + (1. - theta_tumor)*cons_rate_normal)*o2
    solubility = 0.00021 # ml O2 / ml Tissue / kPa
    consfactors = 2000. * 60 * solubility * 1.e2  # 2000. from diff koff, 60 from s to min conversion.
    consumption *= consfactors


    avgStr = lambda t: r'$%s \pm %s$' % tuple(f2l(q) for q in t)

    outfn = 'oxycheck-'+splitext(basename(fn))[0]+'.pdf'
    with mpl_utils.PdfWriter(outfn) as pdfwriter:
      fig = pyplot.figure(figsize=(5,5))

      ax = fig.add_axes([0.05, 0.05, 4.2/5, 4.2/5])
      img = imshow(ax, imslice(o2), ld, vmin = 0., vmax = 9., cmap = matplotlib.cm.Reds, worldscaling = 1.e-3)
      colorbar(ax.get_figure(), ax, img)
      contour(ax, imslice(dist_tumor), ld, levels = [0.], colors = ['k'], worldscaling = 1.e-3)

      tissueavg, tumoravg = tumorAndTissueAverages(o2)
      print 'normal po2: %f +/- %f' %  tissueavg
      print 'tumor po2:  %f +/- %f' %  tumoravg

      fig.suptitle(r'O2 [kPa] (arterial $P_{O_2} = 12 kPa$)')
      fig.text(0.1, 0.90, '<tissue> ' + avgStr(tissueavg) + '\n' +
                          '<tumor> ' + avgStr(tumoravg))

      pdfwriter.savefig(fig)

      fig = pyplot.figure(figsize=(5,5))

      ax = fig.add_axes([0.05, 0.05, 4.2/5, 4.2/5])
      img = imshow(ax, imslice(consumption), ld, vmin = consumption.min(), vmax = consumption.max(), cmap = matplotlib.cm.Reds, worldscaling = 1.e-3)
      colorbar(ax.get_figure(), ax, img)
      contour(ax, imslice(dist_tumor), ld, levels = [0.], colors = ['k'], worldscaling = 1.e-3)

      tissueavg, tumoravg = tumorAndTissueAverages(consumption)
      print 'mro2 normal o2: %f +/- %f' %  tissueavg
      print 'mro2 tumor o2:  %f +/- %f' %  tumoravg

      fig.suptitle(r'$MRO_2$ [ml $O_2$ / 100 ml Tissue min]')
      fig.text(0.1, 0.90, '<tissue> ' + avgStr(tissueavg) + '\n' +
                          '<tumor> ' + avgStr(tumoravg))

      pdfwriter.savefig(fig)

def analyzeO2_without_tumor(f):
    
    dataman = myutils.DataManager(100, [plotBulkTissue2d.DataTumorTissueSingle()])

    
    resultfile = plotBulkTissue2d.ResultFile(f, dataman)
    group = resultfile.groupnames[-1]
    ld = resultfile.obtain_data('ld')

    print group

    arterial_reference = 12.

    theta_tumor = resultfile.obtain_data('fieldvariable', 'theta_tumor', group)
    o2 = arterial_reference * resultfile.obtain_data('fieldvariable','oxy', group)
    dist_tumor = resultfile.obtain_data('fieldvariable', 'dist_tumor', group)

    def tumorAndTissueAverages(field):
      tumor = field[theta_tumor > 0.5]
      normal = field[theta_tumor < 0.5]
      tissueavg = (np.average(normal), np.std(normal))
      tumoravg = (np.average(tumor), np.std(tumor))
      return tissueavg, tumoravg

    cons_rate_normal = float(f['parameters'].attrs['o2_cons_coeff_normal'])
    cons_rate_tumor  = float(f['parameters'].attrs['o2_cons_coeff_tumor'])
    consumption = (theta_tumor * cons_rate_tumor + (1. - theta_tumor)*cons_rate_normal)*o2
    solubility = 0.00021 # ml O2 / ml Tissue / kPa
    consfactors = 2000. * 60 * solubility * 1.e2  # 2000. from diff koff, 60 from s to min conversion.
    consumption *= consfactors


    avgStr = lambda t: r'$%s \pm %s$' % tuple(f2l(q) for q in t)

    outfn = 'oxycheck-'+splitext(basename(fn))[0]+'.pdf'
    with mpl_utils.PdfWriter(outfn) as pdfwriter:
      fig = pyplot.figure(figsize=(5,5))

      ax = fig.add_axes([0.05, 0.05, 4.2/5, 4.2/5])
      img = imshow(ax, imslice(o2), ld, vmin = 0., vmax = 9., cmap = matplotlib.cm.Reds, worldscaling = 1.e-3)
      colorbar(ax.get_figure(), ax, img)
      contour(ax, imslice(dist_tumor), ld, levels = [0.], colors = ['k'], worldscaling = 1.e-3)

      tissueavg, tumoravg = tumorAndTissueAverages(o2)
      print 'normal po2: %f +/- %f' %  tissueavg
      print 'tumor po2:  %f +/- %f' %  tumoravg

      fig.suptitle(r'O2 [kPa] (arterial $P_{O_2} = 12 kPa$)')
      fig.text(0.1, 0.90, '<tissue> ' + avgStr(tissueavg) + '\n' +
                          '<tumor> ' + avgStr(tumoravg))

      pdfwriter.savefig(fig)

      fig = pyplot.figure(figsize=(5,5))

      ax = fig.add_axes([0.05, 0.05, 4.2/5, 4.2/5])
      img = imshow(ax, imslice(consumption), ld, vmin = consumption.min(), vmax = consumption.max(), cmap = matplotlib.cm.Reds, worldscaling = 1.e-3)
      colorbar(ax.get_figure(), ax, img)
      contour(ax, imslice(dist_tumor), ld, levels = [0.], colors = ['k'], worldscaling = 1.e-3)

      tissueavg, tumoravg = tumorAndTissueAverages(consumption)
      print 'mro2 normal o2: %f +/- %f' %  tissueavg
      print 'mro2 tumor o2:  %f +/- %f' %  tumoravg

      fig.suptitle(r'$MRO_2$ [ml $O_2$ / 100 ml Tissue min]')
      fig.text(0.1, 0.90, '<tissue> ' + avgStr(tissueavg) + '\n' +
                          '<tumor> ' + avgStr(tumoravg))

      pdfwriter.savefig(fig)


if __name__ == "__main__":
  f = h5py.File(sys.argv[1], 'r')
  if 'VESSELFILE_MESSAGE' in f.attrs.keys():
    #I think this is a filewithout tumor
    analyzeO2_without_tumor(f)
  else:
    analyzeO2_with_tumor(f)