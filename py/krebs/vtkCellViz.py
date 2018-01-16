#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:39:16 2018

@author: thierry
"""
import vtk

def doit(goodArguments):
  print("loading file: %s" % goodArguments.fileName)
  print("at group: %s" % goodArguments.grpName)
  
  pts = vtk.vtkPoints();
  pts.InsertNextPoint(0,0,0);
  pts.InsertNextPoint(20,0,0);
  pts.InsertNextPoint(60,0,0);
  
  radii = vtk.vtkFloatArray()
  radii.SetName("Radius")
  radii.InsertNextTuple([3.3])
  radii.InsertNextTuple([20.3])
  radii.InsertNextTuple([50.3])
  
  colors = vtk.vtkUnsignedCharArray()
  colors.SetName("Colors")
  colors.SetNumberOfComponents(3)
#
  colors.InsertNextTuple([255, 0, 0])
  colors.InsertNextTuple([0, 255, 0])
  colors.InsertNextTuple([0, 0, 255])
  
  polyData=vtk.vtkPolyData()
  polyData.SetPoints(pts)
  #polyData.GetPointData().SetScalars(colors)
  polyData.GetPointData().SetScalars(radii)
  
  ss = vtk.vtkSphereSource()
  glyph = vtk.vtkGlyph3D()
  glyph.SetScaleModeToScaleByScalar()
  
  glyph.SetSourceConnection(ss.GetOutputPort())
  glyph.SetInputData(polyData)
  glyph.SetInputArrayToProcess(0,0,0,0, 'Radius')
  #print(glyph)
  #glyph.SetInputArrayToProcess(1,0,0,0, 'Colors')
  #glyph.SetColorModeToColorByScalar()
  #glyph.SetSourceConnection(ss.GetOutputPort())
  
  
  
  
  #glyph.ScalingOff()
  #glyph.Update()
  
  #polyData.GetPointData().SetScalars(radii)
  #glyph.SetScaleModeToScaleByScalar()
  #glyph.SetSourceConnection(ss.GetOutputPort())
  glyph.Update()
  
  mapper = vtk.vtkPolyDataMapper()
  mapper.SetInputConnection(glyph.GetOutputPort())
  mapper.SetColorModeToMapScalars()
  mapper.ScalarVisibilityOn()
  mapper.SetScalarRange(0,255)
  mapper.SelectColorArray('Colors')
  
  #mapper.SelectColorArray('Colors')
  # set up the actor
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  
  # do renderer setup stuff
  ren = vtk.vtkRenderer()
  renWin = vtk.vtkRenderWindow()
  renWin.AddRenderer(ren)
  renWin.SetSize(640, 480)
  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)
  
  # add the actor to the renderer
  ren.AddActor(actor)
  
  # render
  iren.Initialize()
  renWin.Render()
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