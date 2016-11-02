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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

from plotDrug import *
import make_video

def generate_frames(filename, outfilename):
  f = h5py.File(filename, 'r+')
  dataman = myutils.DataManager(100, [ DataDrugSingle(), DataDrugAverages(), plotIff.DataTissue(), DataDrugMovie() ])

  times = dataman('movieinfo', f)
  timestr = lambda t: myutils.f2s(t / 3600., latex = True)

  idx_list = range(len(times))
  #idx_list = [0, 1]

  concframes_in = dataman('movieframes', f, idx_list, 'conc_in')
  concframes_in = np.asarray(concframes_in)
  max_conc = np.amax(concframes_in)

  ld = dataman('ld', f)
  tc = dataman('tissue_composition', f)

  filenames = []

  rc = matplotlib.rc
  rc('figure', **{'subplot.left' : 0.01,
                  'subplot.right' : 0.9,
                  'subplot.bottom' : 0.01,
                  'subplot.top' : 0.99,
                  'subplot.wspace' : 0.,
                  'subplot.hspace' : 0.})
  rc('savefig', facecolor = 'white')
  with mpl_utils.SinglePageWriter(outfilename, postfixformat = 'frame%04i') as pdfpages:
    def plt(ax, idx):
      data = concframes_in[idx].transpose()
      #ax.set(title = 't = $%s$' % timestr(times[idx]))
      p1 = imshow(ax, data, ld, vmin=0., vmax=max_conc, cmap=CM.spectral)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      at = mpl_utils.AnchoredText('%s h' % timestr(times[idx]), loc=2, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)


    for idx in idx_list[::2]:
      fig, ax = pyplot.subplots(1, 1, figsize = np.asarray([480, 480]) / 180., dpi = 180)
      print 'frame', idx, 't = ', times[idx]
      mpl_utils.remove_frame(ax)
      plt(ax, idx)
      frame_fn = pdfpages.savefig(fig, dpi = 180)
      pyplot.close()
      filenames.append(frame_fn)
  return filenames

if __name__ == '__main__':
  fn = sys.argv[1]
  outfn = sys.argv[2]
  frame_filenames = generate_frames(fn, outfn)
  make_video.encode(
    pat = outfn+'_frame*.png',
    #files = frame_filenames,
    out = outfn+'.avi', fps = 29.94)