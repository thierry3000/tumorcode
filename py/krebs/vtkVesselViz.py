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

import sys, os
import numpy as np
from vtk import *
from vtkcommon import *
import time

# render vessels with tube filter


def assignedAttribute(dataset, name, attributeType, attributeLoc):
  aa = vtkAssignAttribute()
  if isinstance(dataset, vtkAlgorithmOutput):
    aa.SetInputConnection(dataset)
  else:
    if vtk.VTK_MAJOR_VERSION > 5:
      aa.SetInputData(dataset)
    else:
      aa.SetInput(dataset)
  aa.Assign(name, attributeType, attributeLoc)
  return aa

def getDataSetAttributeBounds(dataset, name, attributeType, attributeLoc):
#  if attributeLoc == 'POINT_DATA':
#    #d = dataset.GetPointData()
#    return dataset.GetScalarRange()
#  elif attributeLoc == 'CELL_DATA':
#    #d = dataset.GetCellData()
#    return dataset.GetCellRange()
  #return d.GetArray(name).GetRange()
  #return d.GetScalars(name).GetRange()
#  abc = dataset.GetCellData()
#  de = abc.GetArray(name)
#  a = dataset.GetCellBounds(attributeType,0)
  return dataset.GetScalarRange()

def getColorTransferFunc(dataspec, bounds):
  name, _, _ = dataspec
  f = vtkColorTransferFunction()
  f.AddRGBPoint( 0, .1, .1, .1)
  f.AddRGBPoint( 1, .9, .9, .9)
  return f

if __name__ == '__main__':
  fn = sys.argv[1]
  fn2 = sys.argv[2]
  ren= vtkRenderer()

  if True:
    reader = vtkPolyDataReader()
    reader.SetFileName(fn)
    reader.Update();
    dataset = reader.GetOutput()
    del reader
    dataset.ComputeBounds()
    vesselbounds = vtkGetDataSetBounds(dataset)

    assign_radius = assignedAttribute(dataset, 'point_radius', 'SCALARS', 'POINT_DATA')

    tube = vtkTubeFilter()
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.SetNumberOfSides(8)
    tube.CappingOff()
    tube.SetGenerateTCoordsToOff()
    tube.SetInputConnection(assign_radius.GetOutputPort())

    dataspec = ('pressure', 'SCALARS', 'POINT_DATA')

    assign_final_attr = assignedAttribute(tube.GetOutputPort(), *dataspec)

    bounds = getDataSetAttributeBounds(dataset, *dataspec)
    bounds = dataset.GetPointData().GetArray(0).GetValueRange()
    mapper = vtkPolyDataMapper()
    mapper.SetScalarRange(*bounds)
    mapper.CreateDefaultLookupTable ()
    mapper.SetInputConnection(assign_final_attr.GetOutputPort())
    mapper.ImmediateModeRenderingOn()

    actor = vtkActor()
    actor.SetMapper(mapper)
    ren.AddActor(actor)

  if True:
    fieldreader = vtkDataSetReader()
    fieldreader.SetFileName(fn2)
    fieldreader.Update()
    fieldataset = fieldreader.GetOutput()
    del fieldreader
    fieldataset.ComputeBounds()

    func = vtkPlane()
    func.SetNormal((0.,0.,1.))
    func.SetOrigin((0.,0.,500))

    cutter = vtkCutter()
    cutter.SetCutFunction(func)
    cutter.SetNumberOfContours(1)
    cutter.SetValue(0, 0.)
    if vtk.VTK_MAJOR_VERSION > 5:
      cutter.SetInputData(fieldataset)
    else:
      cutter.SetInput(fieldataset)
    #cutter.SetInput(fieldataset)

    dataspec = ('fieldOxy', 'SCALARS', 'CELL_DATA')
    bounds = getDataSetAttributeBounds(fieldataset, *dataspec)
    bounds = fieldataset.GetCellData().GetArray(0).GetValueRange()
    print 'fielddata is %s with bounds %s' % (dataspec, bounds)
    
    assign_final_attr = assignedAttribute(cutter.GetOutputPort(), *dataspec)

    colorTransferFunction = getColorTransferFunc(dataspec, bounds)

    mapper = vtkPolyDataMapper()
    mapper.SetScalarRange(*bounds)
    #mapper.CreateDefaultLookupTable ()
    mapper.SetLookupTable(colorTransferFunction)
    mapper.SetInputConnection(assign_final_attr.GetOutputPort())
    mapper.ImmediateModeRenderingOn()

    actor = vtkActor()
    actor.SetMapper(mapper)
    origin = np.min(vesselbounds, axis=0)
    origin[0:2] = 0.
    actor.SetPosition(*origin)
    
    ren.AddActor(actor)
  
  ren.ResetCamera()
  c = ren.GetActiveCamera()
  c.Dolly(1.4)
  c.SetClippingRange(c.GetPosition()[2]*0.9, c.GetPosition()[2]*1.1)

  renWin = vtkRenderWindow()
  renWin.AddRenderer(ren)
  renWin.SetSize(1200, 1200)

  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)

  iren.Start()
  #renWin.Render()

  