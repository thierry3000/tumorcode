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
import posixpath
import numpy as np
import collections
import matplotlib
import matplotlib.pyplot as plt
import extensions
from matplotlib.backends.backend_pdf import PdfPages
import itertools

from mystruct import Struct


def _ncycle(iterable):
    """
    Method to create a generator from an iterable.  It keeps the
    current position of the iterable in memory.  Each time the
    next() method for the iterable is called, it will return the
    next item.  If there are no more items, it will cycle to the
    first item.
    """
  
    for item in itertools.cycle(iterable):
        yield item

markers = _ncycle("osv<>^hdx*")





if __name__ == "__main__":
  fn = sys.argv[1]
  out_fn = os.path.splitext(fn)[0]

  pdf = PdfPages(out_fn + ".pdf")

  f = h5py.File(fn, 'r')
  model_lambda = f.attrs['lambda']
  bytime = collections.defaultdict(list)
  bymethod = collections.defaultdict(list)
  for g in f['/'].itervalues():
    if g.name == '/exact':
      exact = g
    else:
      bytime[g.attrs['dt']].append(g)
      bymethod[g.attrs['method']].append(g)

  for dt, groups in sorted(bytime.iteritems(), key = lambda (k, v): k):
    if dt < exact['x'][-1]*0.01: 
      continue
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.suptitle('dt = %f' % dt)
    plot = fig.add_subplot(1,1,1)
    plot.set_yscale("log")
    plot.set_xlabel('x')
    plot.set_ylabel('y')
    for marker, g in zip(markers, groups):
      plot.plot(np.asarray(g['x']), np.asarray(g['y']), label = g.attrs['method'], marker=marker)
    plot.plot(np.asarray(exact['x']), np.asarray(exact['y']), label = 'exact', linewidth=2)
    plot.legend(loc = 'upper left')
    pdf.savefig(fig)

  err = collections.defaultdict(collections.defaultdict)
  for g in f['/'].itervalues():
      if g.name == '/exact': continue
      val = abs(g['y'][-1] - exact['y'][-1])
      m, t = g.attrs['method'], g.attrs['dt']
      err[m][t] = val
      print m, t, val
  for k, d in err.iteritems():
    err[k] = np.asarray(sorted(d.iteritems(), key = lambda (k, v): k)).transpose()

  fig = plt.figure(figsize=(11.7, 8.3))
  plot = fig.add_subplot(1,1,1)
  plot.set_yscale("log")
  plot.set_xscale("log")
  plot.set_xlabel('dt')
  plot.set_ylabel('e')
  for marker, (method, (t, y)) in zip(markers, err.iteritems()):
    plot.plot(t, y, label = method, marker=marker)
  plot.legend(loc = 'lower right')
  pdf.savefig(fig)

  pdf.close()