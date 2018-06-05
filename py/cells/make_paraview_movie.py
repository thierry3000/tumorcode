# trace generated using paraview version 5.5.0

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def create_screen_shot(number):

  

  fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtk = LegacyVTKReader(FileNames=['/localdisk/thierry/tmp2/fakeTumMTS-default-typeI-sample00-milotti_mts_8_safe_out%04i-cells.vtk'%number])
  
  # get active view
  renderView1 = CreateRenderView('RenderView')
  # uncomment following to set a specific view size
  renderView1.ViewSize = [832, 555]
  
  # show data in view
  fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay = Show(fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtk, renderView1)
  
  # trace defaults for the display properties.
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.Representation = 'Surface'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.ColorArrayName = [None, '']
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.OSPRayScaleArray = 'cell_age'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.SelectOrientationVectors = 'None'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.ScaleFactor = 17.77385025024414
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.SelectScaleArray = 'None'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.GlyphType = 'Arrow'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.GlyphTableIndexArray = 'None'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.GaussianRadius = 0.888692512512207
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.SetScaleArray = ['POINTS', 'cell_age']
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.OpacityArray = ['POINTS', 'cell_age']
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.SelectionCellLabelFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.SelectionPointLabelFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
  
  # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
  fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1399450.0, 1.0, 0.5, 0.0]
  
  # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
  fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1399450.0, 1.0, 0.5, 0.0]
  
  # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.XTitleFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.YTitleFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.ZTitleFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.XLabelFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.YLabelFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.DataAxesGrid.ZLabelFontFile = ''
  
  # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
  #fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
  
  # reset view to fit data
  renderView1.ResetCamera()
  
  #changing interaction mode based on data extents
  renderView1.InteractionMode = '3D'
  
  # create a new 'Glyph'
  glyph1 = Glyph(Input=fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtk,
      GlyphType='Sphere')
  glyph1.Scalars = ['POINTS', 'cell_radii']
  glyph1.Vectors = ['POINTS', 'None']
  glyph1.ScaleMode = 'scalar'
  glyph1.ScaleFactor = 1.0
  glyph1.GlyphMode = 'All Points'
  glyph1.GlyphTransform = 'Transform2'
  
  # Properties modified on glyph1.GlyphType
  glyph1.GlyphType.Radius = 1.0
  
  # show data in view
  glyph1Display = Show(glyph1, renderView1)
  
  # get color transfer function/color map for 'cell_radii'
  cell_ph_exLUT = GetColorTransferFunction('cell_pH_ex')
  cell_ph_exLUT.RescaleTransferFunction(7.3, 7.5)
  
  # trace defaults for the display properties.
  glyph1Display.Representation = 'Surface'
  glyph1Display.ColorArrayName = ['POINTS', 'pH_ex']
  glyph1Display.LookupTable = cell_ph_exLUT
  glyph1Display.OSPRayScaleArray = 'cell_radii'
  glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
  glyph1Display.SelectOrientationVectors = 'None'
  glyph1Display.ScaleFactor = 1.0
  glyph1Display.SelectScaleArray = 'cell_radii'
  #glyph1Display.GlyphType = 'Arrow'
  #glyph1Display.GlyphTableIndexArray = 'cell_radii'
  #glyph1Display.GaussianRadius = 0.9259407043457032
  #glyph1Display.SetScaleArray = ['POINTS', 'cell_radii']
  #glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
  #glyph1Display.OpacityArray = ['POINTS', 'cell_radii']
  #glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
  #glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
  #glyph1Display.SelectionCellLabelFontFile = ''
  #glyph1Display.SelectionPointLabelFontFile = ''
  #glyph1Display.PolarAxes = 'PolarAxesRepresentation'
  
  # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
  glyph1Display.ScaleTransferFunction.Points = [0.9052150249481201, 0.0, 0.5, 0.0, 4.5682501792907715, 1.0, 0.5, 0.0]
  
  # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
  glyph1Display.OpacityTransferFunction.Points = [0.9052150249481201, 0.0, 0.5, 0.0, 4.5682501792907715, 1.0, 0.5, 0.0]
  
  # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
  #glyph1Display.DataAxesGrid.XTitleFontFile = ''
  #glyph1Display.DataAxesGrid.YTitleFontFile = ''
  #glyph1Display.DataAxesGrid.ZTitleFontFile = ''
  #glyph1Display.DataAxesGrid.XLabelFontFile = ''
  #glyph1Display.DataAxesGrid.YLabelFontFile = ''
  #glyph1Display.DataAxesGrid.ZLabelFontFile = ''
  
  # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
  #glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
  #glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
  #glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
  #glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
  
  # show color bar/color legend
  glyph1Display.SetScalarBarVisibility(renderView1, True)
  
  # set scalar coloring
  #ColorBy(glyph1Display, ('POINTS', 'pH_ex'))
  
  # Hide the scalar bar for this color map if no visible data is colored by it.
  HideScalarBarIfNotNeeded(cell_ph_exLUT, renderView1)
  
  # rescale color and/or opacity maps used to include current data range
  #glyph1Display.RescaleTransferFunctionToDataRange(True, False)
  
  # show color bar/color legend
  glyph1Display.SetScalarBarVisibility(renderView1, True)
  
  

  # get opacity transfer function/opacity map for 'pH_ex'
  #pH_exPWF = GetOpacityTransferFunction('pH_ex')
  
  # Rescale transfer function
  #pH_exPWF.RescaleTransferFunction(7.3, 7.5)
  
  # rescale color and/or opacity maps used to exactly fit the current data range
  #glyph1Display.RescaleTransferFunctionToDataRange(False, True)
  
  # rescale color and/or opacity maps used to exactly fit the current data range
  #glyph1Display.RescaleTransferFunctionToDataRange(False, True)
  
  # set scalar coloring
  #ColorBy(glyph1Display, ('POINTS', 'pH_ex'))
  
  
  
  # rescale color and/or opacity maps used to include current data range
  #glyph1Display.RescaleTransferFunctionToDataRange(True, False)
  
  # show color bar/color legend
  #glyph1Display.SetScalarBarVisibility(renderView1, True)
  
  # rescale color and/or opacity maps used to exactly fit the current data range
  #glyph1Display.RescaleTransferFunctionToDataRange(False, True)
  
  # current camera placement for renderView1
  renderView1.CameraPosition = [557.0539332262093, 3.6198501586914062, -0.2259521484375]
  renderView1.CameraFocalPoint = [1.0385475158691406, 3.6198501586914062, -0.2259521484375]
  renderView1.CameraViewUp = [0.0, 0.0, 1.0]
  renderView1.CameraParallelScale = 143.90737119186
  renderView1.Background = [0.0, 0.0, 0.0]
  # Hide orientation axes
  renderView1.OrientationAxesVisibility = 0
  
  # update the view to ensure updated data information
  renderView1.Update()
  
  # save screenshot
  SaveScreenshot('/localdisk/thierry/tmp2/%04i.png' % number, renderView1, ImageResolution=[832, 555])
  
  Delete(glyph1Display)
  del glyph1Display
  Delete(fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtk)
  del fakeTumMTSdefaulttypeIsample00milotti_mts_8_safe_out0455cellsvtk
  # destroy renderView1
  Delete(renderView1)
  del renderView1
if __name__ == '__main__':
  for i in range(001,503):
    print('creating screenshot for timestep %i' %i)
    create_screen_shot(i)