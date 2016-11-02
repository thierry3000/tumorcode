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
from os.path import join, basename, dirname, splitext
import os
import h5py
import h5files
import numpy as np
import collections
import posixpath
import fnmatch
import posixpath
import md5
import cPickle
import uuid
import decimal

import identifycluster

if identifycluster.getname()=='snowden':
  cluster_threads = 16
elif identifycluster.getname()=='durga':
  cluster_threads = 6
else:
  cluster_threads = 4

commonprefix = posixpath.commonprefix
def commonsuffix(l):
  rev = map(lambda s: s[::-1], l)
  ret = commonprefix(rev)
  ret = ret[::-1]
  return ret

def namestr(obj, namespace):
    mybuff = [name for name in namespace if namespace[name] is obj]
    return mybuff[0]

def buildLink(gdst, linkname, gsrc):
  fn = sanitize_posixpath(basename(gsrc.file.filename)) if gsrc else ''
  gdst.attrs[linkname+'_PATH'] = gsrc.name if gsrc else ''
  gdst.attrs[linkname+'_FILE'] = fn
  if gsrc:
    gdst[linkname] = h5py.SoftLink(gsrc.name) if (gsrc.file == gdst.file) else h5py.ExternalLink(fn, gsrc.name)

def splitcommonpresuffix(l):
  '''split common prefix and suffix from a list of strings so that 
     prefix + suffix = input in the edge case of a single input string'''
  prefix  = commonprefix(l)
  restinv = map(lambda s: s[len(prefix):][::-1], l) # remove prefix and reverse char sequence
  suffixinv  = commonprefix(restinv)
  return prefix, suffixinv[::-1]


def updated(d1, d2):
  d3 = type(d1)(d1)
  d3.update(d2)
  return d3


def checksum(*objs):
  '''computes a checksum as md5 of the pickled arguments'''
  return md5.new(cPickle.dumps(objs)).hexdigest()


def uuidstr():
  return str(uuid.uuid1())


def bbitwise_and(flags, bit):
  ''' *b* in bbitwise stands for bool. We want to get an array of type bool
  here which can be used   directly to select a subset of some other array
  (see numpy indexing) '''
  return np.asarray(np.bitwise_and(flags, bit), dtype = np.bool)


def iterate_items(d, keys, skip = False):
  '''iterate over 'key-value pairs for a sequence of given keys'''
  for k in keys:
    try:
      value = d[k]
    except KeyError, e:
      if skip:
        continue
      else:
        raise e
    yield k, value


def strip_from_end(s, end):
  '''removed 'end' from string 's' if s endswith end, otherwise returns s unchanged'''
  return s[:-len(end)] if s.endswith(end) else s


def strip_from_start(s, start):
  '''removed 'start' from string 's' if s startswith start, otherwise returns s unchanged'''
  return s[len(start):] if s.startswith(start) else s


def random_sign(shape, dtype):
  '''returns an array filled with random sign +/-1'''
  return np.asarray(np.sign(np.random.uniform(-1., 1., shape)), dtype = dtype)


def sanitize_posixpath(p):
  return posixpath.normpath(p).strip(posixpath.sep)


def random_stratified(start, end, cnt):
  a = np.linspace(start, end, cnt, endpoint=False)
  r = np.random.uniform(0, (end-start)/cnt, size = len(a))
  return a+r


def processLargeDataset(ds, fun, reduction, blocksize = 1000000, depth=0):
  ''' 
    Circumvent floating point overfloat and accumulation of errors by
    splitting the problem into smaller parts and processing each part
    separately.
    
    Input:
      ds: dataset
      fun: compute some quantity from a subsection of the dataset
      reduction: combine the result of two subsections of the dataset
    Output:
      final call of reduction
  '''
  #print 'pLD'+('  '*depth)+str(ds)
  shape = ds.shape
  if np.product(shape) <= blocksize:  # if the remaining dataset is relatively small
    result = fun(ds)
  else:
    # We split the dataset along the longest axis in half and process each piece separately.
    longest_axis = np.argmax(shape)
    l = shape[longest_axis]
    idx = [slice(None)] * len(shape)
    idx[longest_axis] = slice(0, l/2)
    r1 = processLargeDataset(ds[tuple(idx)], fun, reduction, blocksize, depth+1) # h5py requires a tuple as index list (idx), numpy on the other hand also works with lists
    idx[longest_axis] = slice(l/2, l)
    r2 = processLargeDataset(ds[tuple(idx)], fun, reduction, blocksize, depth+1)
    # 3rd and 4th argument give the fractional volume covered by each of the two sub-boxes
    result = reduction(r1, r2, 1.*(l/2) / l, 1.*(l - l/2) / l)
  return result
  
def largeDatasetAverage(ds, **kwargs):
  def fun(ds):
    return np.average(np.asarray(ds))
  def reduction(a, b, wa, wb):
    return a*wa + b*wb
  return processLargeDataset(ds, fun, reduction, **kwargs)

def largeDatasetAverageAndStd(ds, **kwargs):
  mu = largeDatasetAverage(ds, **kwargs)  
  def fun(ds):
    arr = np.asarray(ds)
    arr = arr - mu
    return np.average(np.square(arr))
  def reduction(a, b, wa, wb):
    return a*wa + b*wb  
  var = processLargeDataset(ds, fun, reduction, **kwargs)
  return mu, np.sqrt(var)



def testProcessLargeDataset():
  def fun(ds):
    print 'fun on ',ds[...]
    return np.average(ds), 1.*np.product(ds.shape)
  
  def reduction(a, b, wa, wb):
    (a,Na),(b,Nb) = a, b
    #wa = Na/(Na+Nb)
    #wb = Nb/(Na+Nb)
    return (a*wa + b*wb), (Na+Nb)

  print 'test1'
  a = np.arange(20)
  result = processLargeDataset(a, fun, reduction, blocksize = 3)
  print result, 'vs', np.average(a)
  
  print 'test2'
  a = np.arange(30).reshape((3,5,2))
  result = processLargeDataset(a, fun, reduction, blocksize = 3)
  print result, 'vs', np.average(a)
  
  print 'test3'
  mu, sd = largeDatasetAverageAndStd(a, blocksize = 3)
  print mu, sd, 'vs', np.average(a), np.std(a)


def closest_items(seq, key, picks):
  """
    for each item in picks find an item
    in seq which is closest to it.
    key is used to get values from items
    in seq for comparison with picks
  """
  from bisect import bisect_left
  l = [ (key(q), q) for q in seq ]
  l.sort(key = lambda q: q[0])
  ltimes = [t for t,q in l]
  result = []
  for t in picks:
    index = bisect_left(ltimes, t)
    assert index>=0
    if index >= len(ltimes):
      result.append(l[-1][1])
    elif index <= 0:
      result.append(l[0][1])
    else:
      # check the object to the left of t
      t0 = ltimes[index]
      t1 = ltimes[index-1]
      if abs(t1-t)<abs(t0-t):
        index-=1
      result.append(l[index][1])
  return result

import pprint as pprint_

class MyPrettyPrinter(pprint_.PrettyPrinter):
  def _format(self, obj, *args, **kwargs):
    if isinstance(obj, collections.defaultdict):
      obj = dict(obj.iteritems())
    elif isinstance(obj, MeanValueArray):
      #obj = (obj.__class__.__name__, np.asarray(obj.avg))
      obj = "%s %s %s" % ( obj.__class__.__name__, obj.sum.dtype, obj.sum.shape)
    return pprint_.PrettyPrinter._format(self, obj, *args, **kwargs)

def pprint(obj):
  p = MyPrettyPrinter()
  p.pprint(obj)


def splitH5PathsFromFilename(fn):
  """
    as a way to specify a path within a h5 file,
    this function requires a string formated like 'filename[:path]'
    and outputs the tuple (filename, path) or (filename, '.')
  """
  i = fn.rfind(':')
  if i > 0:
    (a, b) = fn[:i], fn[i+1:]
    b = b.split(',')
    return (a,b)
  else:
    return fn, ['.']


def splitPath(path):
  def s_(p):
    head, tail = posixpath.split(p)
    if head:
      first = s_(head)
      return first+(tail,)
    else:
      return (tail,)
  if path.startswith(posixpath.sep):    
    return ('',)+s_(path[1:])
  else:
    return s_(path)


class JustHoldsAString(object):
  def __init__(self, t):
    self.t = t
  def __str__(self):
    return self.t


def getTimeSortedGroups(supergroup, namestart = '', key='time'):
  """
    get subgroups sorted by time attribute
    returns list of hdf groups
    startname -> filter out groups the name of which does not start with namestart. Uses the local path name within supergroup. Beware of preceeding '/'
  """
  groups = [ g for k, g in supergroup.iteritems() if ((not namestart or k.startswith(namestart)) and key in g.attrs) ]
  groups.sort(key = lambda g: g.attrs[key])
  return groups



H5FileReference = collections.namedtuple('H5FileReference', ['fn','path'])


def require_snapshot_group_(fmeasure, *args):
  if len(args)==1 and isinstance(args[0], h5py.highlevel.Group):
    h, = args
    return require_snapshot_group_(fmeasure, h.name, h.attrs['time'])
  else:
    name, t = args
    g = fmeasure.require_group(name)
    if not 'time' in g.attrs:
      g.attrs['time'] = t
    return g



def MeasurementFile(f, h5files, prefix='measure_'):
  if not isinstance(f, (str, unicode)):
    fn = f.filename
  else:
    fn = f
  fnpath = dirname(fn)
  fnbase = basename(fn).rsplit('.h5')[0]
  fnmeasure = join(fnpath, prefix+fnbase+'.h5')

  existed = os.path.isfile(fnmeasure)
  fm = h5files.open(fnmeasure, 'a')

  if not existed:
    fm.attrs['TYPE'] = 'MEASUREMENT'
    fm.attrs['SOURCE'] = str(fn)
    fm['source'] = h5py.ExternalLink(fn, '/')
  return fm

# this is how to create a memory file
#      fdst = h5py.File(fnmeasure, 'a', driver = 'core', backing_store = False)


def is_measurement_file(fn):
  with h5py.File(fn, 'r') as f:
    return f.attrs.get('TYPE','') == 'MEASUREMENT'



def hdf_write_dict_hierarchy(grp, name, data):
  """
    store a hierarchy of dicts in a hdf file
    grp: destination hdf group under which to put the data
    name: data name, '.', '', None is also accepted and data is put directly in grp
    data: a dictonary hierarchy. data which is not a dict is copied via grp[name] = data.
  """
  if isinstance(data, dict):
    if name not in ('','.',None):
      g = grp.require_group(name)
    else:
      g = grp
    for k, v in data.iteritems():
      hdf_write_dict_hierarchy(g, k, v)
  elif isinstance(data, np.ma.masked_array):
    ds = grp.create_dataset(name, data=data.torecords())
    ds.attrs['CLASS'] = 'MASKED_ARRAY'
  else:
    #print 'writing %s = %s' % (name, str(data))
    grp[name] = data
    #print 'readback = ', grp[name][...]

def hdf_read_dict_hierarchy(o):
  if isinstance(o, h5py.Group):
    return dict(
      (k, hdf_read_dict_hierarchy(v)) for (k,v) in o.iteritems()
    )
  elif isinstance(o, h5py.Dataset):
    if o.shape == (1,):
      print o.dtype
    a = np.asarray(o)
    if np.ndim(a) == 0:
      a = np.asscalar(a)
    return a

def hdf_read_dict_hierarchy_attr(o):
  assert isinstance(o, h5py.Group)
  d = dict(o.attrs.items())
  for k, v in o.iteritems():
    if isinstance(v, h5py.Group):
      d[k] = hdf_read_dict_hierarchy_attr(v)
  return d


def hdf_data_caching(read, write, f, path, versions = None):
  """
    This function helps with storing and retrieving data from hdf files.

    Example:
      def read_func(g, name):
        return np.asarray(g[name])

      def write_func(g, name):
        g.create_dataset(name, data = some_array)

      hdf_data_caching(read_func, write_func, h5file, ('path','to','data'), (0,0,1))

    If the dataset specified by path exists and have the right versions, then
    read will be called with the last path item as name and the parent hdf group
    as first argument.
    Else the write func will be called first.

    Versioning works by storing an attribute VERSION through the hdf hierarchy.
    It is only stored if the version number is > 0. If the file version is
    smaller than the requested version then the whole group is deleted, the path
    regenerated and write() called to rewrite the data.
  """
  if not isinstance(path, (tuple, list)):
    path = [ path ]
  if versions is None:
    versions = [ 0 ] * len(path)

  def is_version_bad(v, stored_group):
    if v is None: return False # means unversioned and is always good
    try:
      vstored = stored_group.attrs['VERSION']
    except KeyError:
      return True
    return v <> vstored

  def set_metadata(g, v):
    g.attrs['UUID'] = uuidstr()  # unused at the moment. Useful for later
    if v is None: return
    g.attrs['VERSION'] = v

  def output(s):
    #print s
    pass

  p = sanitize_posixpath(path[0])
  v = versions[0]

  if p in f and f[p] == f['.']: # protection against path elements that do nothing
    assert v is None or v == 0
    return hdf_data_caching(read, write, f, path[1:], versions[1:])

  if p in f and is_version_bad(v, f[p]):
    output('* HDF Data Caching - CLEARING %s/%s' % (f.name, p))
    del f[p]

  if len(path)==1:
    if not p in f:
      output('* HDF Data Caching - WRITING %s/%s ...' % (f.name, p))
      write(f, p)
      set_metadata(f[p], v)
      f.file.flush()
      output('* HDF Data Caching - DONE %s/%s ...' % (f.name, p))
    else:
      output('* HDF Data Caching - Reading %s/%s' % (f.name, p))
    ret = read(f, p)
  else:
    if not p in f:
      g = f.create_group(p)
      set_metadata(f[p], v)
    else:
      g = f[p]
    ret = hdf_data_caching(read, write, g, path[1:], versions[1:])

  assert not is_version_bad(v, f[p])
  return ret




def make_hashable(a):
  if isinstance(a, (list, tuple)):
    return tuple(make_hashable(q) for q in a)
  if isinstance(a, dict):
    return tuple((k, make_hashable(q)) for k, q in a.iteritems())
  return a


def UsesDataManager(func):
  '''Wrap function so that it looks up the return value in the cache of
     DataManager before the original function (the argument of UsesDataManager)
     is called'''
  def wrapper(dataman, *args, **kwargs):
    dataman.add_handler_(func) # register using function name
    ret = dataman.obtain_data(func.__name__, *args, **kwargs)  # using the LRU_Cache in dataman
    # will call func if data is not cached
    return ret
  return wrapper
    
    

class DataManager(object):
  def __init__(self, cachesize = 50, handlers = []):
    self.handlers_ = {}    
    for h in handlers:
      self.add_handler_(h)

    def arg2key_(args):
      return make_hashable(args)

    def obtain_data_for_cache_(dataname, *args, **kwargs):
      h = self.handlers_[dataname]
      if DataManager.isDataClass(h):
        return h.obtain_data(self, dataname, *args, **kwargs)
      else:
        return h(self, *args, **kwargs)

    self.cache_ = LRU_Cache(obtain_data_for_cache_, maxsize = cachesize, key = arg2key_)
  #### end of __init__  ###

  def add_handler_(self, h):
    if DataManager.isDataClass(h):
      for k in h.keywords:
        assert not k in self.handlers_
        self.handlers_[k] = h
    else:
      assert (not h.__name__ in self.handlers_) or (self.handlers_[h.__name__] is h)
      self.handlers_[h.__name__] = h

  def __call__(self, *args):
    return self.obtain_data(*args)

  def obtain_data(self, *args, **kwargs):
#    import time
#    t_ = time.time()

    # for testing
    #--------------------
#      dataname, args = args[0], args[1:]
#      h = self.handlers_[dataname]
#      r = h.obtain_data(self, dataname, *args, **kwargs)
    #--------------------
    r = self.cache_(*args, **kwargs)
    #--------------------
#    t_ = (time.time()-t_)*1.e3
#    if t_ > 10.:
#      print '%s, dt = %f ms' % (str(args[0]), t_)
    return r

  @staticmethod
  def isDataClass(h): # as opposed to function
    return hasattr(h, 'keywords') and hasattr(h, 'obtain_data')

#  def __getattr__(self, key):
#    try:
#      return self.__dict__[key]
#    except KeyError:
#      def fwd_(*args, **kwargs):
#        return self.obtain_data(key, *args, **kwargs)
#      return fwd_



# http://stackoverflow.com/questions/4443920/python-building-a-lru-cache
# by
# http://stackoverflow.com/users/1001643/raymond-hettinger
class LRU_Cache(object):
  def __init__(self, original_function, maxsize=1000, key = (lambda x: x)):
    self.original_function = original_function
    self.maxsize = maxsize
    self.mapping = {}
    self.arg2key = key

    PREV, NEXT, KEY, VALUE = 0, 1, 2, 3
    self.head = [None, None, None, None]        # oldest
    self.tail = [self.head, None, None, None]   # newest
    self.head[NEXT] = self.tail

  def __call__(self, *arg):
  #def __call__(original_function, *arg)
    PREV, NEXT, KEY, VALUE = 0, 1, 2, 3
    mapping, head, tail, arg2key = self.mapping, self.head, self.tail, self.arg2key
    sentinel = object()

    key = arg2key(arg)

    link = mapping.get(key, sentinel)
    if link is sentinel:
        #print 'LRU miss %s' % str(key)
        value = self.original_function(*arg)
        if len(mapping) >= self.maxsize:
            oldest = head[NEXT]
            next_oldest = oldest[NEXT]
            head[NEXT] = next_oldest
            next_oldest[PREV] = head
            #print 'LRU discarded %s' % str(oldest[KEY])
            del mapping[oldest[KEY]]
        last = tail[PREV]
        link = [last, tail, key, value]
        mapping[key] = last[NEXT] = tail[PREV] = link
    else:
        #print 'LRU hit %s' % str(key)
        link_prev, link_next, key, value = link
        link_prev[NEXT] = link_next
        link_next[PREV] = link_prev
        last = tail[PREV]
        last[NEXT] = tail[PREV] = link
        link[PREV] = last
        link[NEXT] = tail
    return value


# see
# https://docs.python.org/2/library/decimal.html
def f2s(q, exponential=None, prec=3, latex=False):
  """
    generate nicely formated floating point numbers
    q -> float
    exponential -> True, False, None (means auto)
    prec -> int; the maximal count of significant digits to include
  """
  if isinstance(q, (np.float, np.float32, np.float64)):
    q = float(q)
  if isinstance(q, float): # hack around python2.6 limiation
    try:
      q = q.as_integer_ratio()
    except OverflowError:
      return "inf"
    except ValueError:
      return "NaN"
    q = decimal.Decimal(q[0])/decimal.Decimal(q[1])
  #elif not isinstance(q, decimal.Decimal):
  #  q = decimal.Decimal(float(q))
  else:
    q = decimal.Decimal(q)
  p = q.adjusted() # Used for determining the position of the most significant digit with respect to the decimal point.
  if (abs(p)>3 or exponential==True) and exponential<>False:
    expo = p
    q *= decimal.Decimal(10)**(-expo)
#    print 'expo = ',expo
#    print 'q = ', q
    if q >= decimal.Decimal("9.9999"):
      q /= decimal.Decimal(10)
      expo += 1
    q = f2s(q, exponential=False, prec=prec)

    if latex:
      expo = "10^{%i}" % expo
      if q == "1":
        return expo
      elif q == "-1":
        return "-"+expo
      else:
        return q+"\cdot "+expo
    else:
      expo = "e%i" % expo
      return q + expo
  else:
    q = q.quantize(decimal.Decimal(10)**(p-prec+1))
    qi = q.to_integral()
    if q == qi:
      q = str(qi)
    else:
      q = str(q).rstrip('0')
  return q


def multif2s(*args, **kwargs):
  def turnIntoDecimal(q):
    if isinstance(q, (np.float, np.float32, np.float64)):
      q = float(q)
    if isinstance(q, float): # hack around python2.6 limiation
      try:
        q = q.as_integer_ratio()
      except OverflowError:
        return "inf"
      q = decimal.Decimal(q[0])/decimal.Decimal(q[1])
    else:
      q = decimal.Decimal(q)
    return q

  args = map(turnIntoDecimal, args)
  firstArgExponent = args[0].adjusted()
  prec = kwargs.get('prec', 3)
  
  def formatNumber(q):
    q = q.quantize(decimal.Decimal(10)**(firstArgExponent-prec+1))
    qi = q.to_integral()
    if q == qi:
      q = str(qi)
    else:
      q = str(q)
    #  q = str(q).rstrip('0')
    return q
      
  return map(formatNumber, args)
  

def f2l(*args, **kwargs):
  kwargs['latex'] = True
  return f2s(*args, **kwargs)


#for t in [ 1.e-4, 0.2436e-6, 0.00345, 0.01345, 0.3456, 1.56, 45.631, 303., 543667.6 ]:
#  print "%g -> %s" % (t, f2s(t))


def enumerate_seq(seq):
  """
    Like enumerate but returns a tuple (number, is_first, is_last, item)
    Note:
      In contrast to enumerate() this only takes a sequence as argument because it needs the length. It does not work with iterators.
  """
  t = len(seq)-1
  i = -1
  for q in seq:
    i += 1
    yield (i, i==0, i==t, q)




def zipListOfDicts(s, numpy_output = True):
  """
    Convert a list of dictionaries to a dictionary of lists.
    Can produce numpy arrays instead of lists.
  """
  if not len(s): return {}
  d = dict((k,[]) for k in s[0].iterkeys())
  for q in s:
    for k, v in q.iteritems():
      d[k].append(v)
  if numpy_output:
    for k, v in d.iteritems():
      d[k] = np.asarray(v)
  return d


def walkh5_(gparent, pattern, return_h5objects = False, rec_path_=''):
  head, tail = posixpath.split(pattern)
  #print 'walkh5_(%s) @ %s, %s' % (rec_path_, head, tail)
  if head == '' and tail != '':
    head = tail
    tail = ''
    
  if head == '':
    if return_h5objects:
      #print 'walkh5_ obtain', gparent
      return [gparent]
    else:
      #print 'walkh5_ obtain', rec_path_
      return [rec_path_]
      
  if head == posixpath.sep:
    return walkh5_(gparent.file, tail, return_h5objects, rec_path_ = '/')
  else:
    if not any(k in head for k in '*?'):
      try:
        down = gparent[head]
      except KeyError:
        print 'Error: cannot open path %s/%s' % (gparent.name, head)
        return []
      else:
        return walkh5_(down, tail, return_h5objects, rec_path_ = posixpath.join(rec_path_, head))
    else:
      res = []
      for k in gparent.keys():
        if fnmatch.fnmatch(k, head):
          res += walkh5_(gparent[k], tail, return_h5objects, rec_path_ = posixpath.join(rec_path_, k))
      return res


def walkh5(gparent, pattern, return_h5objects = False, rec_path_=''):
  '''completely new and pretty much untested!!!
     walks a bit faster through h5 trees. supports pattern matching
     of paths in the h5 file, using ?,*. As a little extra we can
     alternative top level paths using |, e.g. "po2/out0000|po2/0006".
     Works only on top level though. ? and * can be used on any level.'''
  patterns = pattern.split('|')
  return sum((walkh5_(gparent, pattern, return_h5objects, rec_path_) for pattern in patterns), [])




def scatter_histogram(xdata, ydata, bins, weights=1., xdata2=None, weights2=None):
  """
    If xdata2 is not supplied:
    computes weighted average and std deviation for each bin:
    i.e. ydata is added to bins defined by 'bins' and the location of
    the data point 'xdata'.

    If xdata2 is supplied:
    The result is the histogram of xdata divided by the histogram of xdata2.
    In both histograms the data points are correspondingly weighted by
    'weights' and 'weights2',
  """
  if isinstance(weights, (float, int)):
    weights = weights * np.ones(xdata.shape, dtype=np.float32)
  if xdata2 is not None:
    if weights2 is None:
      weights2 = 1.
    if isinstance(weights2, (float, int)):
      weights2 = weights2 * np.ones(xdata2.shape, dtype=np.float32)
  else:
    xdata2 = xdata
    weights2 = weights
  val_w, _ = np.histogram(xdata2, bins=bins, weights=weights2)
  val_avg, _ = np.histogram(xdata, bins=bins, weights=weights*ydata)
  val_avg = np.ma.masked_where(val_w<=0., val_avg)
  val_avg /= val_w
  val_sqr, _ = np.histogram(xdata, bins=bins, weights=weights*np.power(ydata, 2.))
  val_sqr = np.ma.masked_where(val_w<=0., val_sqr)
  val_sqr /= val_w
  val_std = np.sqrt(val_sqr - val_avg**2 + 1.e-16)
  return val_avg, val_std, val_sqr


def with_numpy_error_flags(**error_flags):
  def the_decorator(func):
    def wrapper(*args, **kwargs):
      e = np.seterr(**error_flags)
      r = func(*args, **kwargs)
      np.seterr(**e)
      #print 'reset err to %s' % str(e)
      return r
    return wrapper
  return the_decorator


def WeightedAverageStd(a, axis=None, weights=None):
  avg = np.average(a, axis=axis, weights=weights)
  std = np.sqrt(np.average(np.square(a - avg), axis=axis, weights=weights))
  return avg, std


class MeanValueArray(object):
  """
      represents mean and std deviation over an array.


      members:
      cnt : number of additions
      sum : sum of additions
      sqr : sum of sqr of additions

      avg : property returns average
      std : property retuns std deviation

      write (h5grp, name) serialize
      read (h5grp, name) deserialize
      __add__ : add two data obj
      fromHistogram1d()
  """

  def __init__(self, cnt, sum, sqr):
    self.cnt = cnt
    self.sum = sum
    self.sqr = sqr

  @property
  @with_numpy_error_flags(divide = 'ignore', invalid = 'ignore')
  def avg(self):
    return np.ma.masked_invalid(self.sum / self.cnt)

  @property
  @with_numpy_error_flags(divide = 'ignore', invalid = 'ignore')
  def var(self):
    return np.ma.masked_invalid(self.sqr / self.cnt - np.power(self.avg, 2.))

  @property
  @with_numpy_error_flags(divide = 'ignore', invalid = 'ignore')
  def var_mean(self):
    ''' variance of the mean, i.e. estimated deviation of the average from the mean '''
    return np.ma.masked_invalid(self.var / self.cnt)

  @property
  def std(self):
    return np.sqrt(self.var)

  @property
  def std_mean(self):
    ''' standard deviation of the mean, i.e. estimated deviation of the average from the mean '''
    return np.sqrt(self.var_mean)

  def copy(self):
    c = MeanValueArray(None, None, None)
    c += self
    return c

  def __add__(self, other):
    c = MeanValueArray(None, None, None)
    c += self
    c += other
    return c

  def __iadd__(self, other):
    if isinstance(other, MeanValueArray):
      if all(x is None for x in (self.cnt, self.sum, self.sqr)):
        self.cnt = np.copy(other.cnt)
        self.sum = np.copy(other.sum)
        self.sqr = np.copy(other.sqr)
      else:
        self.cnt += other.cnt
        self.sum += other.sum
        self.sqr += other.sqr
      return self
    else:
      s = np.array(other, copy = True)
      if all(x is None for x in (self.cnt, self.sum, self.sqr)):
        self.cnt = np.ones_like(s)
        self.sum = s
        self.sqr = np.square(s)
      else:
        self.cnt += 1
        self.sum += s
        self.sqr += np.square(s)
      return self

  def __imul__(self, other):
    self.sum *= other
    self.sqr *= np.power(other, 2.)
    return self

  def __mul__(self, other):
    c = MeanValueArray(None, None, None)
    c += self
    c *= other
    return c

  def __getitem__(self, key):
    if key in ('cnt', 'sum', 'sqr', 'var', 'std', 'var_mean', 'std_mean'):
      return getattr(self, key)
    if isinstance(key, np.ndarray) and key.dtype == np.bool:
      cnt, sum, sqr = tuple(q[key] for q in [self.cnt, self.sum, self.sqr])
      return MeanValueArray(cnt, sum, sqr)
    if isinstance(key, slice):
      return MeanValueArray(self.cnt[key], self.sum[key], self.sqr[key])
    raise KeyError("dunno how to index %s with %s" % (self.__class__.__name__, str(key)))

  def __len__(self):
    return len(self.sum) if self.sum is not None else 0

  def write(self, h5grp, name, compression = 9):
    g = h5grp.create_group(name)
    g.attrs['TYPE'] = 'BinnedData_V01'
    opt = dict(compression = compression)
    g.create_dataset('cnt', data = self.cnt, **opt)
    g.create_dataset('sum', data = self.sum, **opt)
    g.create_dataset('sqr', data = self.sqr, **opt)
    return g

  @staticmethod
  def read(h5grp, name = None):
    g = h5grp if name is None else h5grp[name]
    d = dict( (q, np.asarray(g[q])) for q in ['cnt','sum','sqr'])
    return MeanValueArray(**d)

  @staticmethod
  def fromHistogram1d(bins, x, y, w = 1.):
    if isinstance(w, (float, int)):
      w = w * np.ones(x.shape, dtype=np.float32)
    cnt, _ = np.histogram(x, bins=bins, weights = w)
    sum, _ = np.histogram(x, bins=bins, weights = w * y)
    sqr, _ = np.histogram(x, bins=bins, weights = w * np.power(y, 2.))
    return MeanValueArray(cnt, sum, sqr)

  @staticmethod
  def inOneBin(y):
    cnt = len(y)
    sum = np.sum(y)
    sqr = np.sum(np.power(y, 2.))
    return MeanValueArray(cnt, sum, sqr)

  @staticmethod
  def empty():
    return MeanValueArray(None, None, None)

  @staticmethod
  def fromArray(a):
    a = np.array(a, copy = True)
    return MeanValueArray(np.ones_like(a), a, np.square(a))

  @staticmethod
  def fromSummation(iterable):
    a = MeanValueArray.empty()
    for b in iterable:
      a += b
    return a


def UpdateHierarchical(d1, d2):
  for k2, v2 in d2.iteritems():
    if isinstance(v2, dict):
      if not k2 in d1:
        d1[k2] = type(v2)()
      UpdateHierarchical(d1[k2], v2)
    else:
      d1[k2] = v2


def TestUpdateHierarchical():
  d1 = dict(a = 5,
            b = 6,
            c = dict(
              d = 1,
              e = 2
            ))
  d2 = dict(a = 100,
            c = dict(
              d = 101,
              f = 102
            ))
  UpdateHierarchical(d1, d2)
  import pprint
  pprint.pprint(d1)



class MultiDict(dict):
  from collections import Iterable
  from collections import defaultdict
  def __init__(self, arg = None):
    self.keysets = MultiDict.defaultdict(set)
    if isinstance(arg, MultiDict.Iterable):
      dict.__init__(self, arg)
      for k, v in dict.iteritems(self):
        for ki in k:
          self.keysets[ki].add(k)
    elif arg is not None and hasattr(arg,'__call__'):
      self.make_default = arg
      dict.__init__(self)
    else:
      dict.__init__(self)

  def __getitem__(self, k):
    assert isinstance(k, tuple)
    return dict.__getitem__(self, k)

  def __missing__(self, k):
    assert isinstance(k, tuple)
    new_item = self.make_default()
    dict.__setitem__(self, k, new_item)
    for ki in k:
      self.keysets[ki].add(k)
    return new_item

  def __setitem__(self, k, v):
    assert isinstance(k, tuple)
    dict.__setitem__(self, k, v)
    for ki in k:
      self.keysets[ki].add(k)

  def iterkeys(self, *k):
    assert isinstance(k, tuple)
    if len(k) == 0: return dict.iterkeys(self)
    sets = tuple(self.keysets[ki] for ki in k)
    intersect = sets[0].intersection(*sets[1:])
    return intersect.__iter__()

  def iteritems(self, *k):
    for ki in self.iterkeys(*k):
      yield (ki,dict.__getitem__(self, ki))

  def items(self, *k):
    return list(self.iteritems(*k))

  def itervalues(self, *k):
    for ki in self.iterkeys(*k):
      yield dict.__getitem__(self, ki)

  def values(self, *k):
    return list(self.itervalues(*k))


def printndarray(a):
  print 'ndim',a.ndim
  print 'shape',a.shape
  print 'size',a.size
  print 'nbytes',a.nbytes
  print 'strides',a.strides
  print 'flags',a.flags
  print 'dtype',a.dtype
  print 'itemsize',a.itemsize


if __name__ == '__main__':
  if 1:
    testProcessLargeDataset()
  if 0:
    TestUpdateHierarchical()
  if 0:
    p = LRU_Cache(ord, maxsize=3)
    for c in 'abcdecaeaa':
        print(c, p(c))
  if 0:
    class Test(object):
      def __init__(self, name):
        self.name = name
        print '%s init' % self.name
      def __enter__(self):
        print '%s enter' % self.name
      def __exit__(self, type, value, traceback):
        print '%s exit' % self.name

    print 'start'
    with RAII(Test('a'), Test('b')) as q:
      print 'acquired before bad'
      print q[0].name
      print q[1].name
      raise RuntimeError('bad')
      print 'fail'
    print 'released'

  if 0:
    f = h5py.File('test.hdf5', driver='core', backing_store=False)
    with f.create_group('gtest') as g:
      g.create_dataset('fubar', data = 1.)
    print g, f

  if 0:
    md = MultiDict(lambda : 'defaultitem!')
    md[1,'shit','t'] = 'test1'
    md[2,'shit','t'] = 'test2'
    md[1,'t'] = 'test3'
    print md.values('t')
    print md.values(1,'t')
    print md.values('shit','t')
    print md['k',]
    md = MultiDict((k, v+'_transformed') for k,v in md.iteritems())
    for k, v in md.iteritems('t',1):
      print k,'=',v
    print '----------'
    for k, v in md.iteritems():
      print k,'=',v