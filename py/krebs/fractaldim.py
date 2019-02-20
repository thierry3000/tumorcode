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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os,sys
from os.path import basename
import time
import h5py
import math
import numpy as np
import time
import posixpath
from mystruct import Struct
import myutils

import extensions
import krebsutils
import qsub


def boxcountImage(arr, boxsize, max_boxsize):
    assert arr.max() == 1 and arr.min()==0
    arr = np.atleast_3d(arr)
    res = []
    while 1:
        # count occupied boxes
        n = np.sum(arr)
        res.append((boxsize,n))

        if 0: # debug out
          print arr.shape, boxsize, n, arr.min(), arr.max()
          q = zplane(arr)
          pyplot.imshow(q, interpolation='nearest', vmin=0, vmax=1)
          pyplot.show()

        if max(arr.shape) <= 10 or boxsize>max_boxsize: break

        boxsize *= 2.

        arr = np.maximum(arr[0::2,:,:], arr[1::2,:,:])
        arr = np.maximum(arr[:,0::2,:], arr[:,1::2,:])
        arr = np.maximum(arr[:,:,0::2], arr[:,:,1::2])

        if arr.shape[0]&1:
          arr = np.append(arr, np.zeros_like(arr[-1:,:,:]), axis=0)
        if arr.shape[1]&1:
          arr = np.append(arr, np.zeros_like(arr[:,-1:,:]), axis=1)
        if arr.shape[2]&1:
          arr = np.append(arr, np.zeros_like(arr[:,:,-1:]), axis=2)

    res = np.asarray(res,dtype=np.double).transpose()
    return res



def zplane(a):
    return np.asarray(a[...,a.shape[2]/2])

#def makeInnerTumorFilter(ds, bordersize, thres):
#    scale = krebs.hdf.get_array3d_scale(ds)
#    dens = krebs.hdf.read_array(ds)
#    thres = thres*reduce(lambda a, b: a*b, scale)/(10.**3)
#    binimage = dens>thres
#    #Image.fromArray(255*zplane(binimage)).save('baseimage.png')
#    filledbg = krebs.flood_fill_3d(binimage)
#    #Image.fromArray(255*zplane(filledbg)).save('filledbg.png')
#    distance = krebs.distancemap(filledbg,scale=scale[0])
#    #Image.fromArray(zplane(distance)).save('distance.png')
#    filter = distance < bordersize
#    binimage[np.nonzero(filter)] = 0
#    filter = np.logical_not(filter)
#    #Image.fromArray(255*zplane(binimage)).save('image.png')
#    return binimage, filter, scale[0]


#def filterTumorVessels(binimage, posMulti, edges, pos):
#    pos = pos/posMulti
#    x, y, z = tuple(pos[...,i] for i in xrange(3))
#    filterPos = binimage[(x,y,z)]  # map from position index to bool
#    #print filterPos
#    # map from edge index to pos index
#    e1 = filterPos[edges[...,0]]
#    e2 = filterPos[edges[...,1]]
#    filter = np.logical_or(e1,e2)
#    #print filter
#    return filter


def makeLd(graph, spacing, relative_offset):
    p = graph.nodes['position']
    pmin = np.amin(p, axis=0)
    pmax = np.amax(p, axis=0)
    s = pmax-pmin
    nc = np.maximum(1, np.asarray(s/spacing, dtype=np.int))
    nc = np.bitwise_and(np.bitwise_not(1), np.bitwise_or(2, nc)) # divisible by 2!
    bb = np.vstack(([0,0,0], nc+1)).T.reshape((6,))
    ld = krebsutils.LatticeDataQuad3d(bb, spacing)
    ld.SetCellCentering([True, True, True])
    ld.SetOriginPosition(pmin + spacing * relative_offset)
    return ld


def vesselsToBinaryImage(graph, ld):
    rho = krebsutils.make_vessel_volume_fraction_field(graph.nodes['position'], graph.edgelist, graph.edges['radius'], ld)
    rho *= ((ld.scale**3) / (5.**2*ld.scale))
    if 0: # debug out
      q = zplane(rho)
      pyplot.imshow(q)
      pyplot.contour(q, [0.5, 1., 1.5])
      pyplot.colorbar()
      pyplot.show()
      #print rho.min(), rho.max()
    np.minimum(rho, 1.00001, rho)
    rho = np.asarray(rho, dtype=np.uint8)
    return rho, ld



def calcVesselBoxCounts(graph, spacing, spacing_max):
  """
    calls c++ function which does boxcounting
    iterates over sequence of increasing box sizes
    boxcounting is based on a rasterized vessel volume distribution
    box is occupied if volume_fraction * volume_factor > volume_threshold
  """
  res = []
  while spacing < spacing_max:
    offsets = [ np.array(tuple(0.5 * (x-.5) for x in index)) for index in np.ndindex((2,2,2)) ]
    cnt = 0
    for offset in offsets:
      ld = makeLd(graph, spacing, offset)
      #print 'boxcounting spacing=%f, gridsize=%i, offset=%s' % (ld.scale, max(ld.GetBox(ravel = True)), str(offset))
      print 'boxcounting spacing=%f, gridsize=%i, offset=%s' % (ld.scale, max(ld.GetBox()), str(offset))
      volume_factor = (ld.scale**3) / (5.**2*ld.scale)
      volume_threshold = 1.
      cnt += krebsutils.calc_vessel_boxcounts(graph.nodes['position'], graph.edgelist, graph.edges['radius'], ld, volume_factor, volume_threshold)
    res.append((ld.scale, cnt/len(offsets)))
    spacing *= math.sqrt(2.)
  return np.asarray(res).transpose()



def calcBoxCounts(data, vesselgroup, dataId, opts):
    """
      boxcounting for various configurations
    """
    if dataId.startswith('tumor'):
      raise RuntimeError('not implemented')
#        ld = krebs.hdf.get_ld(file['lattice'])
#        ds = root['tumor/tc_density']
#        binimage, filter, scale = makeInnerTumorFilter(ds, innerOffset, 0.5)
#        data = boxcountImage(binimage,pixelsize=scale)
    else:
        graph = krebsutils.read_vesselgraph(vesselgroup, ['position', 'flags','radius'])
        flags = graph.edges['flags']
        rad = graph.edges['radius']
        # compute filter
        if dataId=='arteries':
            mask = np.logical_and((flags&krebsutils.ARTERY)>0,rad>10.0)
        elif dataId=='veins':
            mask = np.logical_and((flags&krebsutils.VEIN)>0,rad>10.0)
        elif dataId=='arteriovenous':
            mask = np.logical_and((flags&krebsutils.VEIN)|(flags&krebsutils.ARTERY),rad>10.0)
        elif dataId=='capillaries':
            mask = rad<=6.
        elif dataId=='complete':
            mask = None
        elif dataId=='withintumor':
            mask = flags & krebsutils.WITHIN_TUMOR
        if mask is not None:
          subgraph = graph.get_filtered(edge_indices=np.nonzero(mask))
        else:
          subgraph = graph
        if len(subgraph.edgelist): # if graph is not empty
          print 'computing %s' % dataId
          #krebsutils.set_num_threads(opts.num_threads)
          #do boxcounting
          spacing = opts.spacing
          max_spacing = max(krebsutils.LatticeDataGetWorldSize(makeLd(graph, spacing, 0.)))
          bs, bc = calcVesselBoxCounts(subgraph, spacing, max_spacing)
          data[dataId] = dict(bs = bs, bc = bc)



def boxcountTumorVessels(fn, opts):
  items_to_analyze = ['complete', 'withintumor']

  with h5py.File(fn, 'r') as f:
    groups = myutils.getTimeSortedGroups(f['.'], "out")
    for g in groups[-1:]:
      data = {}
      print '%s : %s' % (fn, g.name)
      for item_to_analyze in items_to_analyze:
        calcBoxCounts(data, g['vessels'], item_to_analyze, opts)

      with myutils.MeasurementFile(fn, h5files, opts.prefix) as fdst:
        dstgroup = myutils.require_snapshot_group_(fdst,g)
        fdst.attrs.create('tumorcode_file_type', data='fractaldim_analyze_bulktum')
        #gdst = fdst.require_snapshot_group_(g)
        #gdst = gdst.recreate_group('fractaldim')
        #myutils.hdf_write_dict_hierarchy(gdst, '.', data)
        myutils.hdf_write_dict_hierarchy(dstgroup, 'fractaldim', data)



def boxcountInitialVessels(fn, opts):
  items_to_analyze = ['arteries', 'veins', 'capillaries', 'complete', 'arteriovenous']

  with h5py.File(fn, 'r') as f:
    data = {}
    for item_to_analyze in items_to_analyze:
      calcBoxCounts(data, f['vessels'], item_to_analyze, opts)

    with myutils.MeasurementFile(fn, h5files, opts.prefix) as fdst:
      gdst = fdst.recreate_group('fractaldim')
      myutils.hdf_write_dict_hierarchy(gdst, '.', data)


class Df:
    def __init__(s, group):
        s.bs = np.asarray(group['bs'])
        s.bc  = np.asarray(group['bc'])
        s.lbs = np.log(s.bs)
        s.lbc  = np.log(s.bc)
        s.params = s.fit()

    def fit(s,sel=None):
        from scipy.optimize import leastsq
        func = lambda p, x, y: (x*p[0]+p[1]-y)
        p, success = leastsq(func, (-2, 0),
                             args=(s.lbs[sel] if sel is not None else s.lbs,
                                   s.lbc[sel]  if sel is not None else s.lbc))
        s.params = p
        return p

    def get(self):
        return -self.params[0]

    def getY(self, params = None):
        if not params:
          params = self.params
        return np.power(self.bs, params[0])*math.exp(params[1])

    def getStr(self):
        params = self.params
        return r'$\propto x^{%s}$' % f2s(params[0])








if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute fractal dimension')  
  parser.add_argument('--plot', action='store_true', default=False)
  #parser.add_argument('--q_local', action='store_true', default=False)
  parser.add_argument('--num_threads', default=1)
  #parser.add_argument('--group', default='vessels')
  parser.add_argument('--prefix', default='fractaldim_')  
  parser.add_argument('--spacing', type=float, default=8.)
  
  parser.add_argument('files', nargs='+', type=str, help='files to calculate')
  parser.add_argument('grp_pattern',help='Where to find the data in the file')  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)  
  
  
#    from optparse import OptionParser
#    parser = OptionParser()
#    parser.add_option('--plot', dest='plot', action='store_true', default=False)
#    parser.add_option('--q-local', dest='q_local', action='store_true', default=False)
#    parser.add_option('--num_threads', dest='num_threads', default=1, type='int')
#    parser.add_option('--group', dest='group', default='vessels')
#    parser.add_option('--prefix', dest='prefix', default='fractaldim_')
#    opts, args = parser.parse_args()
#    opts = Struct( (k,getattr(opts, k)) for k in ['plot', 'q_local', 'num_threads', 'group', 'prefix'])

  if goodArguments.plot == False:
    #goodArgguments.spacing
    #opts['spacing'] = 8.
    if 'vessels' in goodArguments.grp_pattern:
      boxcount_func = boxcountInitialVessels
    elif 'out' in goodArguments.grp_pattern:
      boxcount_func = boxcountTumorVessels

#    for fn in goodArguments.files:
#      if goodArguments.q_local:
#        boxcount_func(fn.name, goodArguments)
#      else:
#        qsub.submit(qsub.Func(boxcount_func, fn.name, goodArguments),
#                    num_cpus = goodArguments.num_threads,
#                    days = 1,
#                    name = 'job_'+basename(fn.name),
#                    change_cwd = True
#                    )
    for fn in goodArguments.files:
      qsub.submit(qsub.Func(boxcount_func, fn, goodArguments),
                  num_cpus = goodArguments.num_threads,
                  days = 1,
                  name = 'job_'+basename(fn),
                  change_cwd = True
                  )

  else: # plotting
    try:
      for fn in goodArguments.files:
        with h5py.File(fn, 'r') as f:
          if not f.attrs.get('tumorcode_file_type') == 'fractaldim_analyze_bulktum':
            raise AssertionError('This in not the right filetype %s!' % fn)
    except Exception, e:
      print e.message
      sys.exit(-1)
      
    dflist = []
    for fn in goodArguments.files:
      with h5py.File(fn, 'r') as f:
        dflist.append(Df(f[goodArguments.grp_pattern + '/fractaldim/withintumor']))

    import matplotlib
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as pyplot
    import mpl_utils

    dpi = 80
    rc = matplotlib.rc
    rc('figure', figsize = list(np.asarray([1024, 768])/dpi), dpi=dpi)
    rc('font', size = 11.)
    rc('axes', titlesize = 13., labelsize = 13.)
    rc('legend', fontsize=13.)
    rc('pdf', compression = 6, fonttype = 42)
    rc('savefig', facecolor='none', edgecolor='none', dpi=dpi)
    rc('font', **{'family':'sans-serif'}) #,'sans-serif':['Helvetica']})
    rc('path', simplify_threshold=0.01)

    df_formula = lambda x : (r'$\propto \epsilon ^ {- %s}$'  % myutils.f2s(x, latex=True))

    fig, ax = pyplot.subplots(1,1)
    for df in dflist:
      ax.plot(df.bs, df.bc, ':kx', label = 'measured' if df is dflist[0] else None)

    avg_y0 = np.average([df.params[1] for df in dflist])
#      ax.plot(dflist[0].bs, dflist[0].getY((-3., avg_y0)), 'r', label = df_formula(3.))
#      ax.plot(dflist[0].bs, dflist[0].getY((-2., avg_y0)), 'g', label = df_formula(2.))
#      ax.plot(dflist[0].bs, dflist[0].getY((-1., avg_y0)), 'b', label = df_formula(1.))
    ax.plot(dflist[0].bs, dflist[0].getY((-3., avg_y0)), '-', color='k', label = '%s, %s, and %s' % (df_formula(1.), df_formula(2.), df_formula(3.)))
    ax.plot(dflist[0].bs, dflist[0].getY((-2., avg_y0)), '-', color='k')
    ax.plot(dflist[0].bs, dflist[0].getY((-1., avg_y0)), '-', color='k')

    df_vals = []
    for df in dflist:
      df.fit(np.logical_and(df.bs > 100, df.bs < 1000.))
      df_vals.append(df.get())
    avg_df = np.average(df_vals)
    std_df = np.std(df_vals)

    avg_y0 = np.average([df.params[1] for df in dflist])

    ax.plot(dflist[0].bs, dflist[0].getY((-avg_df-std_df, avg_y0)), 'r', linewidth = 1., label = (r'$\propto \epsilon ^ {- %s \pm %s }$'  % (myutils.f2s(avg_df, latex=True), myutils.f2s(std_df, latex=True, prec=1)))) #label =  df_formula(avg_df))
    ax.plot(dflist[0].bs, dflist[0].getY((-avg_df+std_df, avg_y0)), 'r', linewidth = 1.)

    ax.set(xscale = 'log', yscale = 'log', xlabel = r'boxsize $\epsilon$', ylabel = r'boxcount', xlim = (10., 10000.))
    ax.grid(linestyle=':', linewidth=0.5, color=(0.7,0.7,0.7))

    ax.legend(loc = 3)

    # dslee gets 2.52 +/- 0.05
    print 'df = %f +/- %f' % (avg_df, std_df)

    pyplot.show()

    fig.subplots_adjust(left = 0.3, right = 1.-0.3, top = 1.-0.3, bottom = 0.3)

    commonpath = os.path.commonprefix(goodArguments.files)
    with mpl_utils.PdfWriter(commonpath+'.pdf') as pdfpages:
      pdfpages.savefig(fig, )


