#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:39:16 2018

@author: thierry
"""
import vtk
import h5py
import numpy as np

def doit(goodArguments):
  print("loading file: %s" % goodArguments.fileName)
  print("at group: %s" % goodArguments.grpName)
  
  f = h5py.File(goodArguments.fileName , 'r')
  g = f[goodArguments.grpName]
  pts = vtk.vtkPoints();
  
  pts.InsertNextPoint(0,0,0);
  pts.InsertNextPoint(20,0,0);
  pts.InsertNextPoint(60,0,0);
  
  radii = vtk.vtkFloatArray()
  radii.SetName("Radius")
  
  radii.InsertNextTuple([3.0])
  radii.InsertNextTuple([3.0])
  radii.InsertNextTuple([3.0])
  
  polyData=vtk.vtkPolyData()
  polyData.SetPoints(pts)
  polyData.GetPointData().AddArray(radii)
  
  #reader = vtkPolyDataReader()
  #reader.SetFileName(fn)
  #reader.Update();
  #dataset = reader.GetOutput()
  #del reader
  #dataset.ComputeBounds()
  #vesselbounds = vtkGetDataSetBounds(dataset)

  aa = vtk.vtkAssignAttribute()
  aa.SetInputData(polyData)
  aa.Assign('Radius', 'SCALARS', 'POINT_DATA')
  
  #assign_radius = vtk.assignedAttribute(polyData, 'Radius', 'SCALARS', 'POINT_DATA')

  tube = vtk.vtkTubeFilter()
  tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
  tube.SetNumberOfSides(8)
  tube.CappingOff()
  tube.SetGenerateTCoordsToOff()
  tube.SetInputConnection(aa.GetOutputPort())

  #dataspec = ('pressure', 'SCALARS', 'POINT_DATA')

  #assign_final_attr = assignedAttribute(tube.GetOutputPort(), *dataspec)
  

  #bounds = getDataSetAttributeBounds(dataset, *dataspec)
  #bounds = dataset.GetPointData().GetArray(0).GetValueRange()
  mapper = vtk.vtkPolyDataMapper()
  #mapper.SetScalarRange(*bounds)
  mapper.CreateDefaultLookupTable ()
  
  #mapper.SetInputConnection(assign_final_attr.GetOutputPort())
  mapper.SetInputConnection(tube.GetOutputPort())
  mapper.ImmediateModeRenderingOn()

  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  ren= vtk.vtkRenderer()
  ren.AddActor(actor)
  ren.ResetCamera()
  c = ren.GetActiveCamera()
  c.Dolly(1.4)
  c.SetClippingRange(c.GetPosition()[2]*0.9, c.GetPosition()[2]*1.1)

  renWin = vtk.vtkRenderWindow()
  renWin.AddRenderer(ren)
  renWin.SetSize(1200, 1200)

  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)

  iren.Start()
  
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='vtk', formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
  #parser.add_argument('AdaptionParamSet')
  #parser.add_argument('--listindex', type=int, help="index of value list" )
  #parser.add_argument('--outputFileFolder', type=str, help="where to store the output, default is working directory")
  parser.add_argument('-f', '--fileName',metavar='fn', type=str,help="vtk file to read")
  parser.add_argument('-g', '--grpName',metavar='fn', type=str,help="vtk file to read")
  #parser.add_argument('--redo', default=False,help="create the optimized network",action="store_true")
  goodArguments, otherArguments = parser.parse_known_args()
  
  doit(goodArguments)