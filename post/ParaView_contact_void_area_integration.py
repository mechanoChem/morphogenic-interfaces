# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
output = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
outputDisplay = Show(output, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'c1'
c1LUT = GetColorTransferFunction('c1')
c1LUT.RGBPoints = [0.00901908, 0.231373, 0.298039, 0.752941, 0.50067104, 0.865003, 0.865003, 0.865003, 0.992323, 0.705882, 0.0156863, 0.14902]
c1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'c1'
c1PWF = GetOpacityTransferFunction('c1')
c1PWF.Points = [0.00901908, 0.0, 0.5, 0.0, 0.992323, 1.0, 0.5, 0.0]
c1PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
outputDisplay.Representation = 'Surface'
outputDisplay.ColorArrayName = ['POINTS', 'c1']
outputDisplay.LookupTable = c1LUT
outputDisplay.SelectTCoordArray = 'None'
outputDisplay.SelectNormalArray = 'None'
outputDisplay.SelectTangentArray = 'None'
outputDisplay.OSPRayScaleArray = 'c1'
outputDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
outputDisplay.SelectOrientationVectors = 'None'
outputDisplay.ScaleFactor = 8.0
outputDisplay.SelectScaleArray = 'c1'
outputDisplay.GlyphType = 'Arrow'
outputDisplay.GlyphTableIndexArray = 'c1'
outputDisplay.GaussianRadius = 0.4
outputDisplay.SetScaleArray = ['POINTS', 'c1']
outputDisplay.ScaleTransferFunction = 'PiecewiseFunction'
outputDisplay.OpacityArray = ['POINTS', 'c1']
outputDisplay.OpacityTransferFunction = 'PiecewiseFunction'
outputDisplay.DataAxesGrid = 'GridAxesRepresentation'
outputDisplay.PolarAxes = 'PolarAxesRepresentation'
outputDisplay.ScalarOpacityFunction = c1PWF
outputDisplay.ScalarOpacityUnitDistance = 4.072584449557349
outputDisplay.OpacityArrayName = ['POINTS', 'c1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
outputDisplay.ScaleTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.992323, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
outputDisplay.OpacityTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.992323, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# hide color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView1, False)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=output)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [40.0, 40.0, 4.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [40.0, 40.0, 4.0]

# set active source
SetActiveSource(slice1)

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'c1']
slice1Display.LookupTable = c1LUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleArray = 'c1'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 8.0
slice1Display.SelectScaleArray = 'c1'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'c1'
slice1Display.GaussianRadius = 0.4
slice1Display.SetScaleArray = ['POINTS', 'c1']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'c1']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0111604, 0.0, 0.5, 0.0, 0.987952, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0111604, 0.0, 0.5, 0.0, 0.987952, 1.0, 0.5, 0.0]

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(output)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# change representation type
outputDisplay.SetRepresentationType('Outline')

# set active source
SetActiveSource(slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 0.01]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# Rescale transfer function
c1LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
c1PWF.RescaleTransferFunction(0.0, 1.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
c1LUT.ApplyPreset('Rainbow Desaturated', True)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'out1'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(c1LUT, renderView1)

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'out1'
out1LUT = GetColorTransferFunction('out1')
out1LUT.RGBPoints = [0.004232961802124083, 0.231373, 0.298039, 0.752941, 0.2679925546335853, 0.865003, 0.865003, 0.865003, 0.5317521474650465, 0.705882, 0.0156863, 0.14902]
out1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'out1'
out1PWF = GetOpacityTransferFunction('out1')
out1PWF.Points = [0.004232961802124083, 0.0, 0.5, 0.0, 0.5317521474650465, 1.0, 0.5, 0.0]
out1PWF.ScalarRangeInitialized = 1

# Rescale transfer function
out1LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
out1PWF.RescaleTransferFunction(0.0, 1.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
out1LUT.ApplyPreset('Rainbow Desaturated', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
out1LUT.ApplyPreset('Rainbow Desaturated', True)

# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on slice1Display
slice1Display.ShowTexturesOnBackface = 0

# Properties modified on slice1Display
slice1Display.ShowTexturesOnBackface = 1

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# reset view to fit data
renderView1.ResetCamera(True)

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=slice1)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# Properties modified on integrateVariables1Display
integrateVariables1Display.Assembly = ''

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# set active view
SetActiveView(renderView1)

# reset view to fit data
renderView1.ResetCamera(True)

# set active view
SetActiveView(spreadSheetView1)

# set active source
SetActiveSource(integrateVariables1)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=integrateVariables1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'contact_void_area_percent'
calculator1.Function = 'c1/(80*80)*100'

# show data in view
calculator1Display = Show(calculator1, spreadSheetView1, 'SpreadSheetRepresentation')

# hide data in view
Hide(integrateVariables1, spreadSheetView1)

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# Properties modified on calculator1Display
calculator1Display.Assembly = ''

SelectIDs(IDs=[-1, 0], FieldType=1, ContainingCells=0)

# clear all selections
ClearSelection()

# set active source
SetActiveSource(calculator1)

SelectIDs(IDs=[-1, 0], FieldType=1, ContainingCells=0)

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2', 'c2']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2', 'c2', 'c1']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2', 'c2', 'c1', 'Points_Magnitude']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2', 'c2', 'c1', 'Points_Magnitude', 'Points']

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Block Number', 'mu1', 'mu2', 'out1', 'out2', 'c2', 'c1', 'Points_Magnitude', 'Points', 'Point ID']

SelectIDs(IDs=[-1, 0], FieldType=1, ContainingCells=0)

# create a new 'Plot Data Over Time'
plotDataOverTime1 = PlotDataOverTime(registrationName='PlotDataOverTime1', Input=calculator1)

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')

# show data in view
plotDataOverTime1Display = Show(plotDataOverTime1, quartileChartView1, 'QuartileChartRepresentation')

# trace defaults for the display properties.
plotDataOverTime1Display.AttributeType = 'Row Data'
plotDataOverTime1Display.UseIndexForXAxis = 0
plotDataOverTime1Display.XArrayName = 'Time'
plotDataOverTime1Display.SeriesVisibility = ['c1 (stats)', 'c2 (stats)', 'contact_void_area_percent (stats)', 'mu1 (stats)', 'mu2 (stats)', 'out1 (stats)', 'out2 (stats)']
plotDataOverTime1Display.SeriesLabel = ['c1 (stats)', 'c1 (stats)', 'c2 (stats)', 'c2 (stats)', 'contact_void_area_percent (stats)', 'contact_void_area_percent (stats)', 'mu1 (stats)', 'mu1 (stats)', 'mu2 (stats)', 'mu2 (stats)', 'out1 (stats)', 'out1 (stats)', 'out2 (stats)', 'out2 (stats)', 'X (stats)', 'X (stats)', 'Y (stats)', 'Y (stats)', 'Z (stats)', 'Z (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotDataOverTime1Display.SeriesColor = ['c1 (stats)', '0', '0', '0', 'c2 (stats)', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'contact_void_area_percent (stats)', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'mu1 (stats)', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'mu2 (stats)', '0.6', '0.3100022888532845', '0.6399938963912413', 'out1 (stats)', '1', '0.5000076295109483', '0', 'out2 (stats)', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'X (stats)', '0', '0', '0', 'Y (stats)', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'Z (stats)', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'N (stats)', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Time (stats)', '0.6', '0.3100022888532845', '0.6399938963912413', 'vtkValidPointMask (stats)', '1', '0.5000076295109483', '0']
plotDataOverTime1Display.SeriesPlotCorner = ['c1 (stats)', '0', 'c2 (stats)', '0', 'contact_void_area_percent (stats)', '0', 'mu1 (stats)', '0', 'mu2 (stats)', '0', 'out1 (stats)', '0', 'out2 (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesLabelPrefix = ''
plotDataOverTime1Display.SeriesLineStyle = ['c1 (stats)', '1', 'c2 (stats)', '1', 'contact_void_area_percent (stats)', '1', 'mu1 (stats)', '1', 'mu2 (stats)', '1', 'out1 (stats)', '1', 'out2 (stats)', '1', 'X (stats)', '1', 'Y (stats)', '1', 'Z (stats)', '1', 'N (stats)', '1', 'Time (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime1Display.SeriesLineThickness = ['c1 (stats)', '2', 'c2 (stats)', '2', 'contact_void_area_percent (stats)', '2', 'mu1 (stats)', '2', 'mu2 (stats)', '2', 'out1 (stats)', '2', 'out2 (stats)', '2', 'X (stats)', '2', 'Y (stats)', '2', 'Z (stats)', '2', 'N (stats)', '2', 'Time (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime1Display.SeriesMarkerStyle = ['c1 (stats)', '0', 'c2 (stats)', '0', 'contact_void_area_percent (stats)', '0', 'mu1 (stats)', '0', 'mu2 (stats)', '0', 'out1 (stats)', '0', 'out2 (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesMarkerSize = ['c1 (stats)', '4', 'c2 (stats)', '4', 'contact_void_area_percent (stats)', '4', 'mu1 (stats)', '4', 'mu2 (stats)', '4', 'out1 (stats)', '4', 'out2 (stats)', '4', 'X (stats)', '4', 'Y (stats)', '4', 'Z (stats)', '4', 'N (stats)', '4', 'Time (stats)', '4', 'vtkValidPointMask (stats)', '4']

# add view to a layout so it's visible in UI
AssignViewToLayout(view=quartileChartView1, layout=layout1, hint=2)

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'c1 (stats)', '0', 'c2 (stats)', '0', 'contact_void_area_percent (stats)', '0', 'mu1 (stats)', '0', 'mu2 (stats)', '0', 'out1 (stats)', '0', 'out2 (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'X (stats)', '1', 'Y (stats)', '1', 'Z (stats)', '1', 'c1 (stats)', '1', 'c2 (stats)', '1', 'contact_void_area_percent (stats)', '1', 'mu1 (stats)', '1', 'mu2 (stats)', '1', 'out1 (stats)', '1', 'out2 (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotDataOverTime1Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'X (stats)', '2', 'Y (stats)', '2', 'Z (stats)', '2', 'c1 (stats)', '2', 'c2 (stats)', '2', 'contact_void_area_percent (stats)', '2', 'mu1 (stats)', '2', 'mu2 (stats)', '2', 'out1 (stats)', '2', 'out2 (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotDataOverTime1Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'c1 (stats)', '0', 'c2 (stats)', '0', 'contact_void_area_percent (stats)', '0', 'mu1 (stats)', '0', 'mu2 (stats)', '0', 'out1 (stats)', '0', 'out2 (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotDataOverTime1Display.SeriesMarkerSize = ['N (stats)', '4', 'Time (stats)', '4', 'X (stats)', '4', 'Y (stats)', '4', 'Z (stats)', '4', 'c1 (stats)', '4', 'c2 (stats)', '4', 'contact_void_area_percent (stats)', '4', 'mu1 (stats)', '4', 'mu2 (stats)', '4', 'out1 (stats)', '4', 'out2 (stats)', '4', 'vtkValidPointMask (stats)', '4']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['c2 (stats)', 'contact_void_area_percent (stats)', 'mu1 (stats)', 'mu2 (stats)', 'out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['contact_void_area_percent (stats)', 'mu1 (stats)', 'mu2 (stats)', 'out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['mu1 (stats)', 'mu2 (stats)', 'out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['mu2 (stats)', 'out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['contact_void_area_percent (stats)', 'out1 (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['contact_void_area_percent (stats)', 'out2 (stats)']

# Properties modified on plotDataOverTime1Display
plotDataOverTime1Display.SeriesVisibility = ['contact_void_area_percent (stats)']

# Properties modified on quartileChartView1
quartileChartView1.LeftAxisTitle = 'Contact void area (percent)'

# Properties modified on quartileChartView1
quartileChartView1.BottomAxisTitle = 'time steps'

# set active view
SetActiveView(renderView1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1480, 1520)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [286.27391759288423, -187.61389367520067, 286.39219628342084]
renderView1.CameraFocalPoint = [40.0, 40.0, 4.0]
renderView1.CameraViewUp = [0.18975182879843167, 0.8386754896499833, 0.5105072639326612]
renderView1.CameraViewAngle = 23.918918918918916
renderView1.CameraParallelScale = 116.48496894979402

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).