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



# remove voids ==================================================

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=output)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'c1']
clip1.Value = 0.50067104

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [40.0, 40.0, 4.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [40.0, 40.0, 4.0]

# Properties modified on clip1
clip1.ClipType = 'Scalar'
clip1.Value = 0.55

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'c1']
clip1Display.LookupTable = c1LUT
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'c1'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 8.0
clip1Display.SelectScaleArray = 'c1'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'c1'
clip1Display.GaussianRadius = 0.4
clip1Display.SetScaleArray = ['POINTS', 'c1']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'c1']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = c1PWF
clip1Display.ScalarOpacityUnitDistance = 3.740628086362419
clip1Display.OpacityArrayName = ['POINTS', 'c1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# hide data in view
Hide(output, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(output)

# set active source
SetActiveSource(output)

# show data in view
outputDisplay = Show(output, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView1, True)

# change representation type
outputDisplay.SetRepresentationType('Outline')

# set active source
SetActiveSource(output)

# set active source
SetActiveSource(clip1)

# rename source object
RenameSource('clip_voids', clip1)

# remove voids ==================================================



# Make 4 views ==================================================

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraFocalDisk = 1.0
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView2, layout=layout1, hint=2)

# set active view
SetActiveView(renderView1)

# split cell
layout1.SplitVertical(1, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.StereoType = 'Crystal Eyes'
renderView3.CameraFocalDisk = 1.0
renderView3.BackEnd = 'OSPRay raycaster'
renderView3.OSPRayMaterialLibrary = materialLibrary1

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView3, layout=layout1, hint=4)


# set active view
SetActiveView(renderView2)

# split cell
layout1.SplitVertical(2, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.AxesGrid = 'GridAxes3DActor'
renderView4.StereoType = 'Crystal Eyes'
renderView4.CameraFocalDisk = 1.0
renderView4.BackEnd = 'OSPRay raycaster'
renderView4.OSPRayMaterialLibrary = materialLibrary1

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView4, layout=layout1, hint=6)

# Make 4 views ==================================================


# View 1 =========================================

# set active view
SetActiveView(renderView1)

# get active source.
clip_voids = GetActiveSource()

# get display properties
clip_voidsDisplay = GetDisplayProperties(clip_voids, view=renderView1)

# get color transfer function/color map for 'c1'
c1LUT = GetColorTransferFunction('c1')
c1LUT.RGBPoints = [0.00901908, 0.231373, 0.298039, 0.752941, 0.50067104, 0.865003, 0.865003, 0.865003, 0.992323, 0.705882, 0.0156863, 0.14902]
c1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'c1'
c1PWF = GetOpacityTransferFunction('c1')
c1PWF.Points = [0.00901908, 0.0, 0.5, 0.0, 0.992323, 1.0, 0.5, 0.0]
c1PWF.ScalarRangeInitialized = 1

# set scalar coloring
ColorBy(clip_voidsDisplay, ('POINTS', 'out1'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(c1LUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip_voidsDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'out1'
out1LUT = GetColorTransferFunction('out1')
out1LUT.RGBPoints = [0.20110436021887063, 0.231373, 0.298039, 0.752941, 0.3680536801094353, 0.865003, 0.865003, 0.865003, 0.535003, 0.705882, 0.0156863, 0.14902]
out1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'out1'
out1PWF = GetOpacityTransferFunction('out1')
out1PWF.Points = [0.20110436021887063, 0.0, 0.5, 0.0, 0.535003, 1.0, 0.5, 0.0]
out1PWF.ScalarRangeInitialized = 1

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView1, False)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
out1LUT.ApplyPreset('Rainbow Desaturated', True)

# Rescale transfer function
out1LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
out1PWF.RescaleTransferFunction(0.0, 1.0)

# find source
output = FindSource('output-*')

# set active source
SetActiveSource(output)

# get display properties
outputDisplay = GetDisplayProperties(output, view=renderView1)

# hide color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView1, False)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# View 1 =========================================

# view 2 ===========================================

# set active view
SetActiveView(renderView2)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get active source.
output = GetActiveSource()

# set active source
SetActiveSource(output)

# show data in view
outputDisplay = Show(output, renderView2, 'UnstructuredGridRepresentation')

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

# show color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera(False)

# change representation type
outputDisplay.SetRepresentationType('Outline')

# hide color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView2, False)

# find source
clip_voids = FindSource('clip_voids')

# set active source
SetActiveSource(clip_voids)

# set active source
SetActiveSource(clip_voids)

# show data in view
clip_voidsDisplay = Show(clip_voids, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip_voidsDisplay.Representation = 'Surface'
clip_voidsDisplay.ColorArrayName = ['POINTS', 'c1']
clip_voidsDisplay.LookupTable = c1LUT
clip_voidsDisplay.SelectTCoordArray = 'None'
clip_voidsDisplay.SelectNormalArray = 'None'
clip_voidsDisplay.SelectTangentArray = 'None'
clip_voidsDisplay.OSPRayScaleArray = 'c1'
clip_voidsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_voidsDisplay.SelectOrientationVectors = 'None'
clip_voidsDisplay.ScaleFactor = 8.0
clip_voidsDisplay.SelectScaleArray = 'c1'
clip_voidsDisplay.GlyphType = 'Arrow'
clip_voidsDisplay.GlyphTableIndexArray = 'c1'
clip_voidsDisplay.GaussianRadius = 0.4
clip_voidsDisplay.SetScaleArray = ['POINTS', 'c1']
clip_voidsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.OpacityArray = ['POINTS', 'c1']
clip_voidsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_voidsDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_voidsDisplay.ScalarOpacityFunction = c1PWF
clip_voidsDisplay.ScalarOpacityUnitDistance = 3.740628086362419
clip_voidsDisplay.OpacityArrayName = ['POINTS', 'c1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_voidsDisplay.ScaleTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_voidsDisplay.OpacityTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView2, False)

# set scalar coloring
ColorBy(clip_voidsDisplay, ('POINTS', 'c2'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(c1LUT, renderView2)

# rescale color and/or opacity maps used to include current data range
clip_voidsDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'c2'
c2LUT = GetColorTransferFunction('c2')
c2LUT.RGBPoints = [0.229825094884741, 0.231373, 0.298039, 0.752941, 0.39531204744237053, 0.865003, 0.865003, 0.865003, 0.560799, 0.705882, 0.0156863, 0.14902]
c2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'c2'
c2PWF = GetOpacityTransferFunction('c2')
c2PWF.Points = [0.229825094884741, 0.0, 0.5, 0.0, 0.560799, 1.0, 0.5, 0.0]
c2PWF.ScalarRangeInitialized = 1

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView2, False)

# Rescale transfer function
c2LUT.RescaleTransferFunction(0.0, 1.0)

# Rescale transfer function
c2PWF.RescaleTransferFunction(0.0, 1.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
c2LUT.ApplyPreset('Rainbow Desaturated', True)

# Hide orientation axes
renderView2.OrientationAxesVisibility = 0

# view 2 ===========================================

# view 3 ===========================================

# set active view
SetActiveView(renderView3)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# find source
output = FindSource('output-*')

# set active source
SetActiveSource(output)

# show data in view
outputDisplay = Show(output, renderView3, 'UnstructuredGridRepresentation')

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

# show color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView3, True)

# reset view to fit data
renderView3.ResetCamera(False)

# change representation type
outputDisplay.SetRepresentationType('Outline')

# hide color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView3, False)

# find source
clip_voids = FindSource('clip_voids')

# set active source
SetActiveSource(clip_voids)

# set active source
SetActiveSource(clip_voids)

# show data in view
clip_voidsDisplay = Show(clip_voids, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip_voidsDisplay.Representation = 'Surface'
clip_voidsDisplay.ColorArrayName = ['POINTS', 'c1']
clip_voidsDisplay.LookupTable = c1LUT
clip_voidsDisplay.SelectTCoordArray = 'None'
clip_voidsDisplay.SelectNormalArray = 'None'
clip_voidsDisplay.SelectTangentArray = 'None'
clip_voidsDisplay.OSPRayScaleArray = 'c1'
clip_voidsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_voidsDisplay.SelectOrientationVectors = 'None'
clip_voidsDisplay.ScaleFactor = 8.0
clip_voidsDisplay.SelectScaleArray = 'c1'
clip_voidsDisplay.GlyphType = 'Arrow'
clip_voidsDisplay.GlyphTableIndexArray = 'c1'
clip_voidsDisplay.GaussianRadius = 0.4
clip_voidsDisplay.SetScaleArray = ['POINTS', 'c1']
clip_voidsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.OpacityArray = ['POINTS', 'c1']
clip_voidsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_voidsDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_voidsDisplay.ScalarOpacityFunction = c1PWF
clip_voidsDisplay.ScalarOpacityUnitDistance = 3.740628086362419
clip_voidsDisplay.OpacityArrayName = ['POINTS', 'c1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_voidsDisplay.ScaleTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_voidsDisplay.OpacityTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView3, False)

# set scalar coloring
ColorBy(clip_voidsDisplay, ('POINTS', 'out1'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(c1LUT, renderView3)

# rescale color and/or opacity maps used to include current data range
clip_voidsDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView3, True)

# get color transfer function/color map for 'out1'
out1LUT = GetColorTransferFunction('out1')
out1LUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.14299999999999996, 0.0, 0.0, 0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.42899999999999994, 0.0, 0.501960784314, 0.0, 0.571, 1.0, 1.0, 0.0, 0.7140000000000001, 1.0, 0.380392156863, 0.0, 0.857, 0.419607843137, 0.0, 0.0, 1.0, 0.878431372549, 0.301960784314, 0.301960784314]
out1LUT.ColorSpace = 'RGB'
out1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'out1'
out1PWF = GetOpacityTransferFunction('out1')
out1PWF.ScalarRangeInitialized = 1

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView3, False)

# Hide orientation axes
renderView3.OrientationAxesVisibility = 0

# flip it ---

# set active source
SetActiveSource(output)

# Properties modified on outputDisplay
outputDisplay.Scale = [1.0, 1.0, -1.0]

# Properties modified on outputDisplay.DataAxesGrid
outputDisplay.DataAxesGrid.Scale = [1.0, 1.0, -1.0]

# Properties modified on outputDisplay.PolarAxes
outputDisplay.PolarAxes.Scale = [1.0, 1.0, -1.0]

# set active source
SetActiveSource(clip_voids)

# Properties modified on clip_voidsDisplay
clip_voidsDisplay.Scale = [1.0, 1.0, -1.0]

# Properties modified on clip_voidsDisplay.DataAxesGrid
clip_voidsDisplay.DataAxesGrid.Scale = [1.0, 1.0, -1.0]

# Properties modified on clip_voidsDisplay.PolarAxes
clip_voidsDisplay.PolarAxes.Scale = [1.0, 1.0, -1.0]

# view 3 ===========================================

# view 4 ===========================================
# set active view
SetActiveView(renderView4)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# find source
output = FindSource('output-*')

# set active source
SetActiveSource(output)

# set active source
SetActiveSource(output)

# show data in view
outputDisplay = Show(output, renderView4, 'UnstructuredGridRepresentation')

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

# show color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView4, True)

# reset view to fit data
renderView4.ResetCamera(False)

# change representation type
outputDisplay.SetRepresentationType('Outline')

# hide color bar/color legend
outputDisplay.SetScalarBarVisibility(renderView4, False)

# find source
clip_voids = FindSource('clip_voids')

# set active source
SetActiveSource(clip_voids)

# set active source
SetActiveSource(clip_voids)

# show data in view
clip_voidsDisplay = Show(clip_voids, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip_voidsDisplay.Representation = 'Surface'
clip_voidsDisplay.ColorArrayName = ['POINTS', 'c1']
clip_voidsDisplay.LookupTable = c1LUT
clip_voidsDisplay.SelectTCoordArray = 'None'
clip_voidsDisplay.SelectNormalArray = 'None'
clip_voidsDisplay.SelectTangentArray = 'None'
clip_voidsDisplay.OSPRayScaleArray = 'c1'
clip_voidsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
clip_voidsDisplay.SelectOrientationVectors = 'None'
clip_voidsDisplay.ScaleFactor = 8.0
clip_voidsDisplay.SelectScaleArray = 'c1'
clip_voidsDisplay.GlyphType = 'Arrow'
clip_voidsDisplay.GlyphTableIndexArray = 'c1'
clip_voidsDisplay.GaussianRadius = 0.4
clip_voidsDisplay.SetScaleArray = ['POINTS', 'c1']
clip_voidsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.OpacityArray = ['POINTS', 'c1']
clip_voidsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
clip_voidsDisplay.DataAxesGrid = 'GridAxesRepresentation'
clip_voidsDisplay.PolarAxes = 'PolarAxesRepresentation'
clip_voidsDisplay.ScalarOpacityFunction = c1PWF
clip_voidsDisplay.ScalarOpacityUnitDistance = 3.740628086362419
clip_voidsDisplay.OpacityArrayName = ['POINTS', 'c1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip_voidsDisplay.ScaleTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip_voidsDisplay.OpacityTransferFunction.Points = [0.00901908, 0.0, 0.5, 0.0, 0.5500000000000002, 1.0, 0.5, 0.0]

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView4, False)

# set scalar coloring
ColorBy(clip_voidsDisplay, ('POINTS', 'c2'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(c1LUT, renderView4)

# rescale color and/or opacity maps used to include current data range
clip_voidsDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView4, True)

# get color transfer function/color map for 'c2'
c2LUT = GetColorTransferFunction('c2')
c2LUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.143, 0.0, 0.0, 0.360784313725, 0.285, 0.0, 1.0, 1.0, 0.429, 0.0, 0.501960784314, 0.0, 0.571, 1.0, 1.0, 0.0, 0.714, 1.0, 0.380392156863, 0.0, 0.857, 0.419607843137, 0.0, 0.0, 1.0, 0.878431372549, 0.301960784314, 0.301960784314]
c2LUT.ColorSpace = 'RGB'
c2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'c2'
c2PWF = GetOpacityTransferFunction('c2')
c2PWF.ScalarRangeInitialized = 1

# hide color bar/color legend
clip_voidsDisplay.SetScalarBarVisibility(renderView4, False)

# Hide orientation axes
renderView4.OrientationAxesVisibility = 0

# flip it ---


# set active source
SetActiveSource(output)

# Properties modified on outputDisplay
outputDisplay.Scale = [1.0, 1.0, -1.0]

# Properties modified on outputDisplay.DataAxesGrid
outputDisplay.DataAxesGrid.Scale = [1.0, 1.0, -1.0]

# Properties modified on outputDisplay.PolarAxes
outputDisplay.PolarAxes.Scale = [1.0, 1.0, -1.0]

# set active source
SetActiveSource(clip_voids)

# Properties modified on clip_voidsDisplay
clip_voidsDisplay.Scale = [1.0, 1.0, -1.0]

# Properties modified on clip_voidsDisplay.DataAxesGrid
clip_voidsDisplay.DataAxesGrid.Scale = [1.0, 1.0, -1.0]

# Properties modified on clip_voidsDisplay.PolarAxes
clip_voidsDisplay.PolarAxes.Scale = [1.0, 1.0, -1.0]

# view 4 ===========================================


#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1606, 1485)

#-----------------------------------
# saving camera placements for views


# current camera placement for renderView1
renderView1.CameraPosition = [264.51448504329136, -170.62528075144297, 266.5706543440295]
renderView1.CameraFocalPoint = [40.0, 40.0, 4.0]
renderView1.CameraViewUp = [0.1897518287984316, 0.8386754896499832, 0.5105072639326612]
renderView1.CameraViewAngle = 13.342318059299192
renderView1.CameraParallelScale = 107.34604859632326

# current camera placement for renderView2
renderView2.CameraPosition = [264.7889324097255, -170.88274990135344, 266.89162171757647]
renderView2.CameraFocalPoint = [40.0, 40.0, 4.0]
renderView2.CameraViewUp = [0.1897518287984316, 0.8386754896499832, 0.5105072639326612]
renderView2.CameraViewAngle = 13.342318059299192
renderView2.CameraParallelScale = 107.47989653721642

# current camera placement for renderView3
renderView3.CameraPosition = [161.58063965397946, -74.05926150431524, 138.18908015380603]
renderView3.CameraFocalPoint = [40.0, 40.0, -4.0]
renderView3.CameraViewUp = [0.1897518287984316, 0.8386754896499832, 0.5105072639326612]
renderView3.CameraViewAngle = 24.932614555256066
renderView3.CameraParallelScale = 56.7097875150313

# current camera placement for renderView4
renderView4.CameraPosition = [161.58063965397946, -74.05926150431524, 138.18908015380603]
renderView4.CameraFocalPoint = [40.0, 40.0, -4.0]
renderView4.CameraViewUp = [0.1897518287984316, 0.8386754896499832, 0.5105072639326612]
renderView4.CameraViewAngle = 24.932614555256066
renderView4.CameraParallelScale = 56.7097875150313

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).


SetActiveView(renderView4)
SetActiveView(renderView3)
SetActiveView(renderView2)
SetActiveView(renderView1)