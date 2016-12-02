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
import numpy as np
import math
import sys
from os.path import join, basename, dirname, splitext
#add path to curvature
sys.path.append(join(dirname(__file__),'..'))
from curvature import estimateCurvatures
import identifycluster
if identifycluster.getname()=='snowden':
  import pyvtk as vtk
else:
  import vtk


def vtkArrayFactory(dtype):
  """return a vtk array type for an element type.
     dtype can a native python or numpy dtype type"""
  if isinstance(dtype, float):
    return vtk.vtkFloatArray
  elif isinstance(dtype, int):
    return vtk.vtkIntArray
  elif isinstance(dtype, str):
    return vtk.vtkStringArray
  elif isinstance(dtype, np.dtype):
    d = dict(i=vtk.vtkIntArray,
             u=vtk.vtkUnsignedIntArray,
             f=vtk.vtkDoubleArray,
             b=vtk.vtkUnsignedCharArray,
             S=vtk.vtkStringArray,
             a=vtk.vtkStringArray,
             U=vtk.vtkStringArray,
             V=vtk.vtkVariantArray)
    try:
      return d[dtype.kind]
    except KeyError:
      pass
  raise RuntimeError("no vtkArray conversion for datatype %s" % dtype)


def asVtkArray(data, name = None, vtk_array_type = None):
  """convert numpy array to vtk array"""
  if vtk_array_type == None:
    vtk_array_type = vtkArrayFactory(data.dtype)
  if len(data.shape) < 2:
    data = np.reshape(data, data.shape+(1,))
  n, nc = data.shape
  a = vtk_array_type()
  a.SetName(name)
  a.SetNumberOfComponents(nc)
  a.SetNumberOfTuples(n)
  if nc > 1:
    for i,q in enumerate(data):
      for j,r in enumerate(q):
        a.SetComponent(i, j, r)
  else:
    for i,q in enumerate(data):
      a.SetValue(i, q)
  return a


def fromVtkArray(a):
  """convert vtk array to numpy array"""
  assert isinstance(a, vtk.vtkAbstractArray)
  typedict = {
    'float' : np.float32,
    'double' : np.float64,
    'int' : np.int32,
    'unsigned char' : np.uint8,
    'unsigned int' : np.uint32,
  }
  s = a.GetDataTypeAsString()
  dtype = typedict[s]
  nc = a.GetNumberOfComponents()
  n = a.GetNumberOfTuples()
  r = np.empty((n, nc) if nc>1 else (n,), dtype=dtype)
  #print r.shape, r.dtype, nc, n, s
#  a.ExportToVoidPointer(r.data)
#  return r
  if nc > 1:
    r = np.zeros((n, nc), dtype=dtype)
    for i in xrange(n):
      for j in xrange(nc):
        r[i,j] = a.GetComponent(i, j)
  else:
    r = np.zeros(n, dtype=dtype)
    for i in xrange(n):
      r[i] = a.GetValue(i)
  return r


def fromVtkPoints(pts):
  a = pts.GetData()
  a = fromVtkArray(a)
  return a



def vtkIter(data, count_member_func_name, get_member_func_name):
  """base iterator for items in a vtk object. It requires the names
     of "GetCount" and "GetElement" member functions"""
  max_index = getattr(data.__class__, count_member_func_name)(data)
  get_func = getattr(data.__class__, get_member_func_name)
  #print max_index
  #print get_func
  index = 0
  while index < max_index:
    yield get_func(data, index)
    index += 1

def vtkIterTuples(data):
  """return Iterator which wraps the GetTupleX members of vtkDataArray"""
  nc = data.GetNumberOfComponents()
  return vtkIter(data, "GetNumberOfTuples", "GetTuple"+str(nc))

def vtkIterArrays(vtkfielddata):
  """return Iterator which wraps the GetArray member of vtkFieldData"""
  return vtkIter(vtkfielddata, "GetNumberOfArrays", "GetArray")

def vtkIterPoints(vtkpoints):
  return vtkIter(vtkpoints, "GetNumberOfPoints", "GetPoint")

def vtkIterCells(vtkcells):
  return vtkIter(vtkcells, "GetNumberOfCells", "GetCell")

def vtkIterCellPointIds(vtkcell):
  return vtkIter(vtkcell, "GetNumberOfPoints", "GetPointId")


def asVtkCellArray(data):
  """convert numpy array to vtkCellArray"""
  data = np.atleast_2d(data)
  n, nc = data.shape
  cells = vtk.vtkCellArray()
  for q in data:
    cells.InsertNextCell(nc)
    for r in q:
      cells.InsertCellPoint(r)
  return cells


def asVtkIdList(data):
  """convert sequence to vtkIdList"""
  r = vtk.vtkIdList()
  r.SetNumberOfIds(len(data))
  for i, d in enumerate(data):
    r.SetId(i, d)
  return r

#def vtkIterIdList(data):
  #return vtkIter(data, "GetNumberOfIds", "GetId")
def fromVtkIdList(idl):
  n = idl.GetNumberOfIds()
  a = np.empty(n, dtype=np.int)
  for i in xrange(n):
    a[i] = idl.GetId(i)
  return a


#def add_vectors(data, name, array_type, arr):
  #data.AddArray(asVtkArray(arr, name, array_type))


#def add_scalars(data, name, array_type, arr):
  #"""This is to add a sequence 'arr' named 'name' to a data array 'data' like polydata.GetPointData(). 'array_type' is the type of the created array, e.g. vtkFloatArray"""
  #data.AddArray(asVtkArray(arr, name, array_type))


def vtkGetDataSetBounds(dataset, mode = '[xyz],[xyz]'):
  """get bounding box as 2x3 numpy array"""
  bounds = np.asarray(dataset.GetBounds())
  bounds = np.reshape(bounds, (3, 2))
  if mode == '[xyz],[xyz]':
    return bounds.transpose()
  elif mode == '[xxyyzz]':
    return bounds.reshape((6,))
  else:
    raise RuntimeError('vtkGetDataSetBounds does not know mode %s' % mode)


def copyPolyDataCells(dataset, cellids):
  """copy cell from the sequence cellids to a new vtkPolyData object"""
  result = vtk.vtkPolyData()
  result.Allocate(dataset, 10000, 1000)
  result.GetCellData().CopyAllocate(dataset.GetCellData(), 0, 0, False)
  result.GetPointData().CopyAllocate(dataset.GetPointData(), 0, 0, False)
  result.CopyCells(dataset, asVtkIdList(cellids), None)
  return result


def vtkImageDataToNumpy(dataset, datanames = None):
  assert isinstance(dataset, vtk.vtkImageData), "dataset must be a vtkImageData"
  arr_res = []
  dim = dataset.GetDataDimension()
  shape_point = np.asarray(dataset.GetDimensions()[:dim])
  assert all(q > 1 for q in shape_point)
  shape_cell = shape_point-1
  celldata = dataset.GetCellData()
  pointdata = dataset.GetPointData()
  arrs = {}
  for q in vtkIterArrays(pointdata):
    arrs[q.GetName()] = (q, pointdata)
  for q in vtkIterArrays(celldata):
    arrs[q.GetName()] = (q, celldata)
  if not datanames:
    datanames = arrs.keys()
  for dataname in datanames:
    a, d = arrs[dataname]
    arr = fromVtkArray(a)
    s = shape_cell if d is celldata else shape_point
    #arr = arr.reshape(s[::-1]).transpose() #WARNING!!!!! i replaced this very important part with the line below!!!
    # aparently the x varies fastest in the memory layout
    arr = arr.reshape(np.concatenate((s,arr.shape[1:])), order='F')
    arr_res.append(arr)
  return arr_res


def vtkImageData(**kwargs):
  origin=np.asarray(kwargs.pop('origin', (0.,0.,0.)))
  spacing=np.asarray(kwargs.pop('spacing', (1.,1.,1.)))
  shape = np.ones((3,), dtype=np.int)
  if 'point_shape' in kwargs:
    s = np.asarray(kwargs.pop('point_shape'))
    dim = np.argwhere(s > 1).max()+1
  elif 'cell_shape':
    s = kwargs.pop('cell_shape')
    dim = np.argwhere(s > 1).max()+1
    s = np.asarray(s)+1
  shape[:dim] = s[:dim]
  origin[dim:] = 0
  assert 0 < dim <= 3
  grid = vtk.vtkImageData()
  grid.SetDimensions(shape)
  grid.SetSpacing(*spacing)
  grid.SetOrigin(*origin);
  assert grid.GetDataDimension() == dim
  return grid


def vtkImageDataFromLd(ld, scalefactor=1.):
  import krebsutils
  if type(ld) == krebsutils.LatticeDataQuad3d:
    o = scalefactor * np.asarray(ld.worldBox[[0, 2, 4]])
    h = scalefactor * np.asarray((ld.scale,)*3)
    s = ld.shape
  else:
    #depricated
#    s = np.asarray((ld['SIZEX'], ld['SIZEY'], ld['SIZEZ']))
#    h = scalefactor * np.asarray((ld['SCALE'],)*3)
#    o = scalefactor*(np.asarray(ld['WORLD_OFFSET']))
#    cc = ld['CENTERING']
#    o -= cc * h * 0.5
    s = np.asarray(ld.shape)
    o = np.asarray(ld.GetOriginPosition())
    h = scalefactor * np.asarray((ld.GetScale(),)*3)
    #h = 30.
  return vtkImageData(origin=o, spacing=h, cell_shape=s)


def vtkImageDataAddData(grid, data, at, name):
  dim = grid.GetDataDimension()
  d = grid.GetCellData() if at=='CellData' else grid.GetPointData()
  data = np.atleast_3d(data)  
  if len(data.shape)==4 and data.shape[3] == dim:
    data = np.asarray(data).transpose().reshape((dim,-1)).transpose()
    d.AddArray(asVtkArray(data, name))
  elif len(data.shape) == 3:
    data = np.asarray(data).transpose().flatten()
    d.AddArray(asVtkArray(data, name))
  else:  
    raise RuntimeError("bad field shape: %s" % str(data.shape))


def vtkComputeCellVolume(c):
  dim = c.GetCellDimension()
  idlist = vtk.vtkIdList()
  pts = vtk.vtkPoints()
  c.Triangulate(0, idlist, pts)
  idlist = fromVtkIdList(idlist)
  idlist = idlist.reshape((-1, dim+1))
  pts = fromVtkPoints(pts)
  pts = pts.reshape((-1, dim+1, 3))
  #print idlist
  res = 0.
  if dim == 2:
    for a, b, c in pts:
      vol = vtk.vtkTriangle.TriangleArea(a, b, c)
      res += vol
  else:
    raise RuntimeError("not implemented")
  return res


def vtkIntegrateData(dataset):
  dim = dataset.GetCell(0).GetCellDimension()
  id_buffer = vtk.vtkIdList()
  pts_buffer = vtk.vtkPoints()
  weight_sum = 0.
  
  class ValueHandler:
    def __init__(self, fd):
      self.fd = fd
      self.fd_arrays = [ fd.GetArray(i) for i in xrange(fd.GetNumberOfArrays()) ]
      self.value = np.zeros(sum(a.GetNumberOfComponents() for a in self.fd_arrays), dtype = np.double)
    def addinto(self, i, weight):
      v = self.value
      k = 0
      for a in self.fd_arrays:
	for c in xrange(a.GetNumberOfComponents()):
	  v[k] += weight * a.GetComponent(i, c)
	  k += 1
    def get_values(self):
      values = []
      k = 0
      for a in self.fd_arrays:
	nc = a.GetNumberOfComponents()
	if nc == 1:
	  values.append(self.value[k])
	  k += 1
	else:
	  r = np.empty(nc, dtype = self.value.dtype)
	  for i in xrange(nc):
	    r[i] = self.value[k]
	    k += 1
	  values.append(r)
      return values
      
  pd_handler = ValueHandler(dataset.GetPointData())
  cd_handler = ValueHandler(dataset.GetCellData())
  
  #integrator = vtk.vtkCellIntegrator()  
  
  for cell_index in xrange(dataset.GetNumberOfCells()):
    c = dataset.GetCell(cell_index)
    #weight = integrator.integrate(dataset, cell_index)
    ok = c.Triangulate(0, id_buffer, pts_buffer)  
    if not ok: continue
    idlist = fromVtkIdList(id_buffer)
    idlist = idlist.reshape((-1, dim+1))
    pts = fromVtkPoints(pts_buffer)
    pts = pts.reshape((-1, dim+1, 3))
    weight = 0.
    if dim == 2:
      for x, y, z in pts:
        w = vtk.vtkTriangle.TriangleArea(x, y, z)
        weight += abs(w)
#    elif dim == 1:
#      for x, y in pts:
#        weight += math.sqrt(np.sum(np.square(y-x)))
    else:
      raise RuntimeError("not implemented dim %i" % dim)
    weight_sum += weight
    pd_weight = weight/c.GetNumberOfPoints()
    for i in vtkIterCellPointIds(c):
      pd_handler.addinto(i, pd_weight)
    cd_handler.addinto(cell_index, weight)

  return cd_handler.get_values(), pd_handler.get_values(), weight_sum
  #print "cd:", cd_handler.get_values()
  #print "pd:", pd_handler.get_values()


def fromVtkTriangleMesh(vtkds):
  """
    returns (points, triangles) numpy arrays
  """
  assert isinstance(vtkds, vtk.vtkPolyData)
  triangles = []
  for cell in vtkIterCells(vtkds):
    if cell.GetCellType() <> vtk.VTK_TRIANGLE:
      raise RuntimeError("non triangle cell encountered")
    t = [ i for i in vtkIterCellPointIds(cell) ]
    triangles.append(t)
  points = [ vtkds.GetPoint(i) for i in xrange(vtkds.GetNumberOfPoints()) ]
  triangles = np.asarray(triangles)
  points = np.asarray(points)
  return points, triangles


def vtkGetLineComponents(surf):
  """vtkds must be polydata describing lines
     returns list of index arrays, point coordinates"""
  assert isinstance(surf, vtk.vtkPolyData)
  components = []
  num_points = surf.GetNumberOfPoints()
  mask = np.zeros(num_points, dtype=np.int)
  component_num = 0
  for main_p in xrange(num_points):
    if mask[main_p] <> 0: continue
    l = []
    last_p = -1
    p = main_p
    component_num += 1
    #while mask[p] == 0:
    while True:
      l.append(p)
      if mask[p]<>0:
        break
      mask[p] = component_num
      pts = []
      r = vtk.vtkIdList()
      surf.GetPointCells(p, r)
      #print "at %i discovered: " % p,
      for ci in range(r.GetNumberOfIds()):
        c = surf.GetCell(r.GetId(ci))
        for i in vtkIterCellPointIds(c):
          #print "%i(%s)" % (i, ci),
          if not (i==p or i==last_p):
            pts.append(i)
      #for ci, c in [ (r.GetId(i), surf.GetCell(r.GetId(i))) for i in range(r.GetNumberOfIds()) ]:
        #for i in tuple(c.GetPointId(j) for j in range(2)):
      #print ""
      if not pts:
        break
      last_p = p
      p = pts[0]
    if l:
      components.append(l)
  #print mask
  #pprint.pprint(components)
  points = np.asarray([surf.GetPoint(j)[:2] for j in xrange(num_points)])
  components = [ np.asarray(c, dtype=np.int) for c in components ]
  return components, points


def vtkCurvature(dataset_in, num_points):
  dataset = vtk.vtkPolyData()
  dataset.ShallowCopy(dataset_in)
  components, points = vtkGetLineComponents(dataset)
  curv = estimateCurvatures(components, points, num_points)
  arr = asVtkArray(curv, "curvature", vtk.vtkDoubleArray)
  dataset.GetPointData().AddArray(arr)
  return dataset


def vtkContour(dataset, values):
  """
    computes a contour lines (surfaces) using the vtkContourFilter
    and returns a vtkPolyData object
  """
  #dataset.GetPointData().SetActiveScalars(name)
  iso = vtk.vtkContourFilter()
  iso.SetInputData(dataset)
  if not hasattr(values, '__iter__'):
    values = [values]
  for i, v in enumerate(values):
    iso.SetValue(i, v)
  iso.Update()
  iso = iso.GetOutput()
#  cl = vtk.vtkCleanPolyData()
#  cl.SetInput(iso.GetOutput())
#  cl.SetTolerance(0.00001)
#  cl.Update()
#  iso = cl.GetOutput()
  return iso


def pyContourPolyLines(contour):
  #vtkds.GetPointData().SetActiveScalars(dataname)
#  iso = vtkContour(vtkds, value)
  components, points = vtkGetLineComponents(contour)
  for i, l in enumerate(components):
    components[i] = points[l, :]#np.asarray(points[j] for j in l)
  return components


def vtkCellDataToPointData(dataset):
  #f = vtk.vtkPointDataToCellData()
  f = vtk.vtkCellDataToPointData()
  f.PassCellDataOff()
  f.SetInputData(dataset)
  f.Update()
  return f.GetOutput()
  

def ZippedOpen(fn):
  import gzip
  if fn.endswith('.gz'):
    return gzip.open(fn, 'r')
  else:
    return open(fn,'r')


def ZippedRead(Reader, fn):
  import gzip
  import os
  import tempfile
  org_fn = fn
  delete_fn = False
  if fn.endswith('.gz'):
    org_fn = os.path.splitext(fn)[0]
    f = gzip.open(fn, 'r')
    fd, fn2 = tempfile.mkstemp(suffix='.vtu', dir=os.path.dirname(fn), text=True)
    print 'unzipping to %s' % fn2
    os.fdopen(fd, 'w').writelines(f)
    fn = fn2
    delete_fn = True
  org_fn = os.path.splitext(os.path.basename(org_fn))[0]
  try:
    r = Reader()
    r.SetFileName(fn)
    r.Update()
    if r.GetErrorCode() <> 0:
      raise RuntimeError("%s failed" % str(Reader))
  finally:
    if delete_fn:
      os.remove(fn)
  return r.GetOutput(), org_fn


def vtkDatasetFromHdf5(g):
  """
    extract a vtkDataset Object from a text file representation stored in a
    hdf5 character array dataset g. The concrete returned type depends on the
    content of the file.
  """
  gridtype2vtureader = {
    'vtkImageData' : vtk.vtkXMLImageDataReader,
    'vtkPolyData' : vtk.vtkXMLPolyDataReader,
    'vtkUnstructuredGrid' : vtk.vtkXMLUnstructuredGridReader,
  }
  filetype = g.attrs["TYPE"]
  if filetype == 'VTK_FILE': # shit! i used vtu files at first before i noticed that the legacy vtk format is much easier to deal with ...
    reader = vtk.vtkDataSetReader()
    reader.ReadFromInputStringOn()
    reader.SetInputString(np.asarray(g).tostring())
    reader.Update()  
  else:
    probable_vtk_dataset_type = g.attrs.get('VTK_DATASET_TYPE', 'vtkUnstructuredGrid')
    import mkstemp  # my stuff
    tf = mkstemp.File(suffix=".vtx")
    f = open(tf.filename, 'wb')
    f.write(np.asarray(g).tostring())
    f.close()
    reader = gridtype2vtureader[probable_vtk_dataset_type]()
    reader.SetFileName(tf.filename)
    reader.Update()
    tf.remove()    
  return reader.GetOutput()
  


def vtkDatasetToHdf5(g, name, ds):
  """
    store a vtk dataset ds, formated as vtk file, within a hdf character array dataset.
    The hdf dataset is created with name name under the group g.
  """
  #assert False # lol this is just copy pasted from some old code. needs to be done properly
  #assert type(ds) is vtk.vtkPolyData 
  assert ds.GetClassName() == 'vtkPolyData' # there needs to be a map from dataset type to writer type. Also the python type does not equal the wrapped c++ type.
  writer = vtk.vtkPolyDataWriter()
  writer.SetInput(ds)
  writer.WriteToOutputStringOn()
  writer.SetFileTypeToASCII()
  writer.Update()
  iso = writer.GetOutputString()  
  l = writer.GetOutputStringLength()
  iso = iso[:l] # iso is a string, but only the content up to l is valid!
  ds = g.create_dataset(name, data = np.fromstring(iso, dtype=np.uint8), compression = 'gzip', compression_opts = 9)
  ds.attrs["TYPE"] = "VTK_FILE"
  ds.attrs["VTK_DATASET_TYPE"] = type(ds).__name__



def vtkPlane(normal=(0.,0.,1.), origin=(0.,0.,0.)):
  p = vtk.vtkPlane()
  p.SetNormal(*normal)
  p.SetOrigin(*origin)
  return p
   

def vtkCutDataSet(dataset, func, values):
  """
    uses vtkCutter,
    always creates a mesh consisting only of triangles (in3d)
  """
  f = vtk.vtkCutter()
  f.SetCutFunction(func)
  if not hasattr(values, '__iter__'):
    values = [values]
  for i, v in enumerate(values):
    f.SetValue(i, v)
  f.SetInput(dataset)
  f.Update()
  return f.GetOutput()


def vtkScaleDataSet(dataset, scalefactor):
  if dataset.IsA('vtkImageData'):
    res = vtk.vtkImageData()
    res.ShallowCopy(dataset)
    p = np.asarray(res.GetOrigin())
    p *= scalefactor
    res.SetOrigin(*p)
    p = np.asarray(res.GetSpacing())
    p *= scalefactor
    res.SetSpacing(*p)
  elif dataset.IsA('vtkPointSet'):
    f = vtk.vtkTransformFilter()
    t = vtk.vtkTransform()
    t.Scale(scalefactor,scalefactor,scalefactor)
    f.SetTransform(t)
    f.SetInput(dataset)
    f.Update()
    res = f.GetOutput()
  else:
    raise RuntimeError('not implemented')
  return res
  
  
def matplotCmToLt(cm):
  """turn a matplotlib colormap into a vtkLookupTable"""
  N = cm.N
  c = vtk.vtkLookupTable()
  c.Allocate(N, N)
  for i in xrange(N):
    c.SetTableValue(i, *cm(i))
  return c


def npImageLayout(img):
  """
    my conversion routines between numpy and vtkImageData are consistent, meaning that
    numpy -> vtk -> numpy is the identity operation. BUT my numpy data layout is different
    from the standard layout for bitmap images. This functions converts my layout to 
    the standard layout. 
    It flips the y axis and swaps the row and column axes.
  """
  return img[:,::-1,...].swapaxes(0,1)


def vtkRender2d(dataset, cmap, crange, pixelsize, return_numpy = False):
  """
    uses rendering
  """
  bounds = vtkGetDataSetBounds(dataset, mode='[xyz],[xyz]')  
  size = np.max(np.abs(bounds), axis=0)
  
  res = np.asarray(size*2./pixelsize, dtype=int)
  parallel_scale = size[1]
  
  ren= vtk.vtkRenderer()

  lt = matplotCmToLt(cmap)
  
  rmin, rmax = crange

  mapper = vtk.vtkDataSetMapper()
  mapper.SetInput(dataset)
  mapper.InterpolateScalarsBeforeMappingOn()
  mapper.UseLookupTableScalarRangeOff()
  mapper.SetScalarRange(rmin, rmax)
  mapper.SetLookupTable(lt)

  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  #actor.GetProperty().SetRepresentationToWireframe()
  ren.AddActor(actor)  

  c = ren.GetActiveCamera()
  c.ParallelProjectionOn()
  c.SetParallelScale(parallel_scale)
  c.SetClippingRange(c.GetPosition()[2]*0.9, c.GetPosition()[2]*1.1)

  renWin = vtk.vtkRenderWindow()
  renWin.AddRenderer(ren)
  renWin.SetSize(res[0], res[1])

  c = vtk.vtkWindowToImageFilter()
  c.SetInput(renWin)
  c.Update()

  if return_numpy:
    ds = c.GetOutput()
    img, = vtkImageDataToNumpy(ds)
    return npImageLayout(img)
  else:
    return c.GetOutput()
  
#def resample(vtkds):
#  loc = vtk.vtkCellLocator()
#  loc.SetDataSet(vtkds)
#  loc.Update()
#  print loc.FindClosestPoint((0.,0.,0.) ) # NOT IMPLEMENTED!!