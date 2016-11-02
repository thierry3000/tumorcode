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
import os
import numpy as np

import matplotlib

import identifycluster
if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
  matplotlib.use('agg')

#from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, AuxTransformBox, VPacker,\
     TextArea, DrawingArea

class PdfWriter(object):
  def __init__(self, filename, multipage=True, absorbExtension_ = True):
    import matplotlib.backends.backend_pdf as be
    self.be = be
    if filename.endswith('.pdf') and absorbExtension_:
      self.filename = filename = os.path.splitext(filename)[0]
    if multipage:
      self.pdf = be.PdfPages(filename+'.pdf')
    else:
      self.pdf = None

  def savefig(self, fig, postfix='', **kwargs):
    """
      postfix is added to the filename if multipage is False, else everything is in one file and postfix is ignored
      other mpl options:
        dpi,
        facecolor,
        edgecolor
        orientation -> 'portrait', 'landscape',
        papertype,
        format,
        transparent,
        ...
        see figure.savefig()
    """
    c = self.be.FigureCanvasPdf(fig)
    if self.pdf:
      self.pdf.savefig(fig, **kwargs)
    else:
      c.print_figure(self.filename+postfix+'.pdf', **kwargs)

  def close(self):
    if self.pdf:
      self.pdf.close()

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()




class SinglePageWriter(object):
  def __init__(self, filename, fileformat, postfixformat = None, absorbExtension_ = True):
    self.filename = filename
    self.pagenum = 0
    self.postfixformat = postfixformat if postfixformat else 'page%02i'
    self.fileFormat = fileformat if fileformat is not None else 'svg'
    if filename.endswith('.'+fileformat) and absorbExtension_:
      self.filename = filename = os.path.splitext(filename)[0]

  def savefig(self, fig, postfix = '', **kwargs):
    self.pagenum += 1
    basename, _ = os.path.splitext(self.filename)
    postfix = postfix if postfix else self.postfixformat % self.pagenum
    filename = '%s_%s.%s' % (basename, postfix, self.fileFormat)
    fig.savefig(filename, **kwargs)

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    pass

# note: strange bug, svg has to come before pdf in the list of file formats, or else some exception is thrown
class PageWriter(object):
  def __init__(self, filename, postfixformat = None, fileformats = ['pdf']):
    if any(filename.endswith('.'+fmt) for fmt in fileformats):
      filename = os.path.splitext(filename)[0]
    def WriterFactory(fmt):
      if fmt == 'pdf':
        return PdfWriter(filename, absorbExtension_ = False)
      else:
        return SinglePageWriter(filename, fmt, postfixformat, absorbExtension_ = False)
    self.writers = map(WriterFactory, fileformats)
  
  def savefig(self, *args, **kwargs):
    for w in self.writers:
      w.savefig(*args, **kwargs)
  
  def __enter__(self):
    for w in self.writers:
      w.__enter__()
    return self
  
  def __exit__(self, type, value, traceback):
    for w in self.writers:
      w.__exit__(type, value, traceback)


markers = '.,ov^<>12348sp*hH+XDd|_'

def styleMarker(marker, color = 'k', linewidth = 1.0):
  ''' 
    returns a dict with arguments for the plot function which result in nicely styled (color wise) markers.
    Produces markers with colored outline and transparent interior. More options may be added later.
  '''
  if marker in 'x+':
    return dict(marker = marker, markeredgewidth = linewidth, color = color, markeredgecolor = color)
  else:
    return dict(marker = marker, markeredgewidth = linewidth, color = color, markeredgecolor = color, markerfacecolor = 'none')


def styleMarkerSP(marker, color = 'k', linewidth = 1.0):
  ''' 
    for scatter plots ...
  '''
  if marker in 'x+':
    return dict(marker = marker, linewidth = linewidth, color = color, edgecolor = color)
  else:
    return dict(marker = marker, linewidth = linewidth, color = color, edgecolor = color, facecolor = 'none')




def add_contourlines(axes, lines, **kwargs):
  """
    lines - sequence of sequences of 2d points
    kwargs are passed to add_collection
  """
  return axes.add_collection(matplotlib.collections.LineCollection(lines, transOffset=axes.transData, **kwargs))


def add_contourlines_vtk(axes, dataset, **kwargs):
  """
    dataset - a polydata with line segments
    kwargs are passed to add_collection
  """
  import vtkcommon
  lines = vtkcommon.pyContourPolyLines(dataset)
  return add_contourlines(axes, lines, **kwargs)


#def imshow_vtk(axes, vtkds, cellsize, cmap, crange=None, cscale=1.):
#  import vtkcommon
#  bounds_xxyy = vtkcommon.vtkGetDataSetBounds(vtkds, mode='[xxyyzz]')[:4]
#  data = vtkcommon.fromVtkArray(vtkds.GetCellData().GetScalars())
#  if crange is None:
#    crange = (data.min(), data.max())
#  elif crange == 'zero-centered':
#    q = np.abs(data).max()
#    crange = (-q,q)
#  img = vtkcommon.vtkRender2d(vtkds, cmap, crange, cellsize, return_numpy = True)
#  return axes.imshow(img, extent=bounds_xxyy, interpolation='nearest', cmap=cmap, vmin=cscale*crange[0], vmax=cscale*crange[1])

def MakeTextPage(lines, figsize=None):
  fig, ax = matplotlib.pyplot.subplots(1,1, figsize=figsize)
  text = '\n'.join(lines)
  ax.text(0.01, 0.95, text, verticalalignment = 'top')
  remove_frame(ax)
  ax.axis('off')
  return fig, ax



def log_bins(a, b, g=2.):
  import math
  import sys
  assert a>0 and b>a and b < sys.float_info.max
  assert g>1.
  x = 1.
  imin = 0
  while x > a:
    x = math.pow(g, imin)
    imin -= 1
  imax = imin
  while x < b:
    x = math.pow(g, imax)
    imax += 1
  r = np.arange(imin, imax-imin, 1)
  r = np.power(g, r)
  assert r[0]<a and r[-1]>b
  return r


lc = 'bgrcmy'
ls = ('-','--','-.',':')

def default_lc(i):
  return lc[i%len(lc)]

def default_ls(i):
  return ls[i%len(ls)]


def tight_layout(fig, **kwargs):
  '''
    lets you use tight_layout on a specific fig. Will change the current figure
    ,call pyplot.tight_layout, and change the figure back again.
  '''
  backup = matplotlib.pyplot.gcf()
  matplotlib.pyplot.figure(fig.number)
  matplotlib.pyplot.tight_layout(**kwargs)
  matplotlib.pyplot.figure(backup.number)


def hist(plt, x, **kwargs):
  """
    matplotlib hist plot with custom parameters and scalar weight parameter
  """
  w = kwargs.pop('weights', None)
  if w is not None:
    if not isinstance(w, (list, tuple, np.ndarray)):
      w = w * np.ones(x.shape)
    kwargs['weights'] = w
  if not 'normed' in kwargs:
    kwargs['normed'] = False
  if not 'histtype' in kwargs:
    kwargs['histtype'] = 'step'
#  if 'selection' in kwargs:
#    indx = kwargs['selection']
#    x = x[indx]
#    if w is not None:
#      w = w[indx]
#      kwargs['weights'] = w
  return plt.hist(x, **kwargs)



def errorbar(ax, x, y, **kwargs):
  """
    Wrapper around errorbar to plot white-filled markers.
    Can be used for line plots, too for consistency
  """
  color = kwargs.get('color', None)
  try:
    every = kwargs.pop('every')
  except KeyError:
    pass
  else:
    kwargs['markevery'] = every
    kwargs['errorevery'] = every
  return ax.errorbar(x, y,
              markeredgecolor = kwargs.pop('markeredgecolor', color),
              markerfacecolor = kwargs.pop('markerfacecolor', 'white'),
              capsize = kwargs.pop('capsize', 0),
              **kwargs)

def subplots_adjust_abs(fig, **kw):
  '''adjust plot margins in inches'''
  unit = kw.pop('unit', None)
  unit2factor = { 'mm'   : mm_to_inch,
                  'inch' : None,
                  None   : None}
  factor = unit2factor[unit]
  w,h = fig.get_figwidth(), fig.get_figheight()
  for k in 'left bottom right top wspace hspace'.split():
    try:
      q = kw.pop(k)
    except KeyError:
      continue
    if factor is not None:
      q *= factor
    if k in 'left right wspace'.split():
      q /= w
    else:
      q /= h
    if q < 0.:
      q = 1.+q
    kw[k] = q
    #print k,q
  fig.subplots_adjust(**kw)


def subplots_adjust(fig, **kw):
  abs = kw.pop('abs', False)
  if abs:
    subplots_adjust_abs(fig, **kw)
  else:
    fig.subplots_adjust(**kw)


def add_crosshair(ax, (x, y), **kwargs):
  '''Adds a crosshair through the coordinate point (x,y). kwargs are
     passed to ax.plot. Attempts to keep the axes ranges, and plots
     only as far as the coordinate axes go.'''
  x0, x1 = ax.get_xlim()
  y0, y1 = ax.get_ylim()
  ll = []
  if x is not None:
    x0 = min(x0, x)
    x1 = max(x1, x)
    ll += ax.plot([x, x], [y0, y1], 'k', **kwargs)
  if y is not None:
    y0 = min(y0, y)
    y1 = max(y1, y)
    ll += ax.plot([x0, x1], [y, y], 'k', **kwargs)
  if ll:
    ax.set(xlim = (x0,x1), ylim=(y0,y1))
    for l in ll:
      l.set_zorder(1)
  return ax


def remove_frame(*axes):
  '''removes the frame of a plot including the axes. Useful for
     plotting just some text on a blank page'''
  for ax in axes:
    ax.xaxis.set(visible=False)
    ax.yaxis.set(visible=False)
    ax.set(xticks = [], yticks = [])


mm_to_inch = 0.03937
a4size = np.asarray((210 * mm_to_inch, 297 * mm_to_inch))
a4size_mm = np.asarray((210, 297))

class loc(object):
  '''location constants for matplotlib'''
  upper_right = 1
  upper_left  = 2
  lower_left = 3
  lower_right = 4
  right = 5
  center_left = 6
  center_right = 7
  lower_center = 8
  upper_center = 9
  center = 10


def AccumulateXYBounds(xmin, xmax, ymin, ymax, previous = None):
  if previous is None:
    return (xmin, xmax), (ymin, ymax)
  else:
    (xpmin, xpmax), (ypmin, ypmax) = previous
    return (min(xmin, xpmin), max(xmax, xpmax)), (min(ymin, ypmin), max(ymax, ypmax))



def SetSensibleAxLimits(ax, xbounds, ybounds):
  arg = {}
  for i, (a,b) in enumerate([xbounds, ybounds]):
    m = b-a
    if m > 0: # might be zero which is when we let mpl handle this
      d = 0.05 * m
      a -= d
      b += d
      arg[['xlim','ylim'][i]] = (a,b)
  if arg: # if there is anything in there
    ax.set(**arg)


def add_sizebar(ax, size = 200, text = ur"200 \u03BCm", color = 'white', loc = 4, myfontsize='normal', size_vertical=1.0):
  '''adds a size bar into an image plot'''
  a = matplotlib.font_manager.FontProperties(size=myfontsize)
  asb = AnchoredSizeBar(ax.transData, size, text, loc=4, pad=0.1, borderpad = 0.5, sep=5, frameon=False,fontproperties=a,size_vertical=size_vertical)
  asb.get_child().get_children()[1].get_children()[0].set(color = color)
  asb.get_child().get_children()[0].get_children()[0].set(color = color)
  ax.add_artist(asb)



def subplots_abs_mm(figsize, nrows, ncols, bottom, left, width, height, wspace, hspace, **subplot_args):
  '''like pyplot.subplots but takes margins in mm as arguments'''
  size = np.asarray(figsize) * mm_to_inch
  relmm = lambda s, x: (x * mm_to_inch)/s
  width = relmm(size[0], width)
  height = relmm(size[1], height)
  left = relmm(size[0], left)
  bottom = relmm(size[1], bottom)
  wspace = relmm(size[0], wspace)
  hspace = relmm(size[1], hspace)
  fig = matplotlib.pyplot.figure(figsize = size)
  axes = np.zeros((nrows, ncols), dtype = np.object)
  offseti = height + hspace
  offsetj = width + wspace
  for (i,j) in np.ndindex((nrows, ncols)):
    axes[i,j] = fig.add_axes([offsetj*j+left, offseti*(nrows-i-1)+bottom, width, height], **subplot_args)
  return fig, axes




def testplotexactimage_():
  import numpy as np
  import matplotlib.pyplot as pyplot
  import matplotlib
  matplotlib.use('Cairo')
  """
    plot an image as image file such that the image is mapped identically
  """
  dpi = 70.
  res = 400
  mpl_res = res-1 #int(res*1.293)-2
  im = np.zeros((res, res), dtype=np.uint8)
  for (x,y) in np.ndindex(res,res):
    if x&1 and y&1:
      im[(x,y)] = 1

  fig = pyplot.figure(figsize=(float(mpl_res)/dpi, float(mpl_res)/dpi), dpi=dpi)
  plot = fig.add_subplot(1,1,1)
  plot.yaxis.set_visible(False)
  plot.xaxis.set_visible(False)
  plot.set_axis_off()
  img = plot.imshow(im)
  img.set_interpolation('nearest')
  fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
  fig.savefig("testplotexactimage.png", dpi=dpi) #, bbox_inches='tight', pad_inches=0)



def testVtkRenderPixelExact_():
  import vtkcommon
  import vtk
  import matplotlib.pyplot as pyplot

  res = (10, 20, 0)
  im = np.zeros(res[:2], dtype=np.uint8)
  for (x,y) in np.ndindex(*res[:2]):
    if x&1 and y&1:
      im[(x,y)] = 1
  for (x,y) in np.ndindex(res[0], 1):
    im[(x,res[1]-1)] = 1

  ds = vtkcommon.vtkImageData(cell_shape=res, spacing = (30., 30., 30.), origin=-np.asarray(res)*30.*0.5)
  vtkcommon.vtkImageDataAddData(ds, im, "CellData", "testdata")
  ds.GetCellData().SetActiveScalars("testdata")
  x0, x1, y0, y1, _, _ = vtkcommon.vtkGetDataSetBounds(ds, mode='[xxyyzz]')

  back_im, = vtkcommon.vtkImageDataToNumpy(ds)

  print 'numpy -> vtk -> numpy is identity operation:', np.all(im == back_im)

  writer = vtk.vtkDataSetWriter()
  writer.SetInput(ds)
  writer.SetFileName("plottest-org.vtk")
  writer.Update()

  vtkimg = vtkcommon.vtkRender2d(ds, matplotlib.cm.spectral, (-1, 2), 30.)

  writer = vtk.vtkDataSetWriter()
  writer.SetInput(vtkimg)
  writer.SetFileName("plottest-img.vtk")
  writer.Update()

  writer = vtk.vtkPNGWriter()
  writer.SetInput(vtkimg)
  writer.SetFileName("plottest-img.png")
  writer.Write()

  img, = vtkcommon.vtkImageDataToNumpy(vtkimg)

  pdfpages = PdfPages("plottest-img.pdf")
  fig = pyplot.figure()
  pyplot.subplot(121)
  pyplot.imshow(vtkcommon.npImageLayout(img), extent = (x0, x1, y0, y1), interpolation='nearest')
  pyplot.subplot(122)
  pyplot.imshow(vtkcommon.npImageLayout(im), extent = (x0, x1, y0, y1), interpolation='nearest')
  pdfpages.savefig(fig)
  pdfpages.close()





def testSVGRasterization_():
  #rc = matplotlib.rc
  #rc('svg', image_inline = False)
  import matplotlib.pyplot as pyplot

  img = np.asarray([[1, 2], [3, 4]])

  x = np.random.normal(0., 1., 1000)
  y = np.random.normal(0., 1., 1000)

  fig, axes = pyplot.subplots(1, 2, figsize = (8, 4))

  axes[0].imshow(img)
  axes[1].scatter(x, y, s=1., rasterized=True, edgecolors='none')

  fig.savefig('image_and_raster_test_300dpi.svg', dpi = 300)
  fig.savefig('image_and_raster_test_72dpi.svg', dpi = 72)



def testSVGGreek_():
  import matplotlib.pyplot as pyplot

  fig = pyplot.figure(figsize = (4, 4))

  fig.text(0.1, 0.5, ur'Normal \u03BC $\mu$ \u0040', size = 20)

  fig.savefig('textmu.svg')


def testErrorbar_():
  import matplotlib.pyplot as pyplot

  a = np.linspace(0., 1., 50)
  b = a*a
  c = 0.21 * np.ones_like(b) * np.average(b)

  fig, ax = pyplot.subplots(1,1, figsize = (2, 2), dpi = 300)
  plotit(a, b, fmt = '-s', label = 'fubar', color = 'r', every = 10)
  plotit(a, b+ 3*c, yerr = c, fmt = '-s', label = 'fubar', color = 'm', every = 5)
  plotit(a, b+ 6*c, fmt = ':x', label = 'test3', color = 'g', every = 10)
  ax.legend()

  fig.savefig('test.svg')
  fig.savefig('test.pdf')
  pyplot.show()


if __name__ == '__main__':
  testErrorbar_()

