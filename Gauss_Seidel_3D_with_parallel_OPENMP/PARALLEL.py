# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
gS_3D_WITH_PARALLEL_BLACK_REDcsv = CSVReader(FileName=['/home/mayank/HPC/ASSIGNMENT1/Gauss_Seidel_3D_with_parallel_OPENMP/GS_3D_WITH_PARALLEL_BLACK_RED.csv'])
gS_3D_WITH_PARALLEL_BLACK_REDcsv.DetectNumericColumns = 1
gS_3D_WITH_PARALLEL_BLACK_REDcsv.UseStringDelimiter = 1
gS_3D_WITH_PARALLEL_BLACK_REDcsv.HaveHeaders = 1
gS_3D_WITH_PARALLEL_BLACK_REDcsv.FieldDelimiterCharacters = ','
gS_3D_WITH_PARALLEL_BLACK_REDcsv.AddTabFieldDelimiter = 0
gS_3D_WITH_PARALLEL_BLACK_REDcsv.MergeConsecutiveDelimiters = 0

# Properties modified on gS_3D_WITH_PARALLEL_BLACK_REDcsv
gS_3D_WITH_PARALLEL_BLACK_REDcsv.HaveHeaders = 0

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.UseCache = 0
spreadSheetView1.ViewSize = [400, 400]
spreadSheetView1.CellFontSize = 9
spreadSheetView1.HeaderFontSize = 9
spreadSheetView1.SelectionOnly = 0
spreadSheetView1.GenerateCellConnectivity = 0
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.SelectedComponent = -1
spreadSheetView1.InvertOrder = 0
spreadSheetView1.BlockSize = 1024L
spreadSheetView1.HiddenColumnLabels = []
spreadSheetView1.FieldAssociation = 'Point Data'

# get layout
#layout1 = GetLayout()

# place view in the layout
#layout1.AssignView(2, spreadSheetView1)

# show data in view
gS_3D_WITH_PARALLEL_BLACK_REDcsvDisplay = Show(gS_3D_WITH_PARALLEL_BLACK_REDcsv, spreadSheetView1)

# trace defaults for the display properties.
gS_3D_WITH_PARALLEL_BLACK_REDcsvDisplay.CompositeDataSetIndex = [0]

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Table To Structured Grid'
tableToStructuredGrid1 = TableToStructuredGrid(Input=gS_3D_WITH_PARALLEL_BLACK_REDcsv)
tableToStructuredGrid1.WholeExtent = [0, 0, 0, 0, 0, 0]
tableToStructuredGrid1.XColumn = 'Field 0'
tableToStructuredGrid1.YColumn = 'Field 0'
tableToStructuredGrid1.ZColumn = 'Field 0'

# Properties modified on tableToStructuredGrid1
tableToStructuredGrid1.WholeExtent = [0, 99, 0, 99, 0, 99]
tableToStructuredGrid1.YColumn = 'Field 1'
tableToStructuredGrid1.ZColumn = 'Field 2'

# show data in view
tableToStructuredGrid1Display = Show(tableToStructuredGrid1, spreadSheetView1)

# trace defaults for the display properties.
tableToStructuredGrid1Display.CompositeDataSetIndex = [0]

# hide data in view
Hide(gS_3D_WITH_PARALLEL_BLACK_REDcsv, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1060, 538]

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToStructuredGrid1)

# show data in view
tableToStructuredGrid1Display_1 = Show(tableToStructuredGrid1, renderView1)

# trace defaults for the display properties.
tableToStructuredGrid1Display_1.Representation = 'Outline'
tableToStructuredGrid1Display_1.AmbientColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.ColorArrayName = [None, '']
tableToStructuredGrid1Display_1.DiffuseColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.LookupTable = None
tableToStructuredGrid1Display_1.MapScalars = 1
tableToStructuredGrid1Display_1.MultiComponentsMapping = 0
tableToStructuredGrid1Display_1.InterpolateScalarsBeforeMapping = 1
tableToStructuredGrid1Display_1.Opacity = 1.0
tableToStructuredGrid1Display_1.PointSize = 2.0
tableToStructuredGrid1Display_1.LineWidth = 1.0
tableToStructuredGrid1Display_1.RenderLinesAsTubes = 0
tableToStructuredGrid1Display_1.RenderPointsAsSpheres = 0
tableToStructuredGrid1Display_1.Interpolation = 'Gouraud'
tableToStructuredGrid1Display_1.Specular = 0.0
tableToStructuredGrid1Display_1.SpecularColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.SpecularPower = 100.0
tableToStructuredGrid1Display_1.Luminosity = 0.0
tableToStructuredGrid1Display_1.Ambient = 0.0
tableToStructuredGrid1Display_1.Diffuse = 1.0
tableToStructuredGrid1Display_1.EdgeColor = [0.0, 0.0, 0.5]
tableToStructuredGrid1Display_1.BackfaceRepresentation = 'Follow Frontface'
tableToStructuredGrid1Display_1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.BackfaceOpacity = 1.0
tableToStructuredGrid1Display_1.Position = [0.0, 0.0, 0.0]
tableToStructuredGrid1Display_1.Scale = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.Orientation = [0.0, 0.0, 0.0]
tableToStructuredGrid1Display_1.Origin = [0.0, 0.0, 0.0]
tableToStructuredGrid1Display_1.Pickable = 1
tableToStructuredGrid1Display_1.Texture = None
tableToStructuredGrid1Display_1.Triangulate = 0
tableToStructuredGrid1Display_1.UseShaderReplacements = 0
tableToStructuredGrid1Display_1.ShaderReplacements = ''
tableToStructuredGrid1Display_1.NonlinearSubdivisionLevel = 1
tableToStructuredGrid1Display_1.UseDataPartitions = 0
tableToStructuredGrid1Display_1.OSPRayUseScaleArray = 0
tableToStructuredGrid1Display_1.OSPRayScaleArray = 'Field 3'
tableToStructuredGrid1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tableToStructuredGrid1Display_1.OSPRayMaterial = 'None'
tableToStructuredGrid1Display_1.Orient = 0
tableToStructuredGrid1Display_1.OrientationMode = 'Direction'
tableToStructuredGrid1Display_1.SelectOrientationVectors = 'Field 3'
tableToStructuredGrid1Display_1.Scaling = 0
tableToStructuredGrid1Display_1.ScaleMode = 'No Data Scaling Off'
tableToStructuredGrid1Display_1.ScaleFactor = 0.192156
tableToStructuredGrid1Display_1.SelectScaleArray = 'Field 3'
tableToStructuredGrid1Display_1.GlyphType = 'Arrow'
tableToStructuredGrid1Display_1.UseGlyphTable = 0
tableToStructuredGrid1Display_1.GlyphTableIndexArray = 'Field 3'
tableToStructuredGrid1Display_1.UseCompositeGlyphTable = 0
tableToStructuredGrid1Display_1.UseGlyphCullingAndLOD = 0
tableToStructuredGrid1Display_1.LODValues = []
tableToStructuredGrid1Display_1.ColorByLODIndex = 0
tableToStructuredGrid1Display_1.GaussianRadius = 0.0096078
tableToStructuredGrid1Display_1.ShaderPreset = 'Sphere'
tableToStructuredGrid1Display_1.CustomTriangleScale = 3
tableToStructuredGrid1Display_1.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
tableToStructuredGrid1Display_1.Emissive = 0
tableToStructuredGrid1Display_1.ScaleByArray = 0
tableToStructuredGrid1Display_1.SetScaleArray = ['POINTS', 'Field 3']
tableToStructuredGrid1Display_1.ScaleArrayComponent = ''
tableToStructuredGrid1Display_1.UseScaleFunction = 1
tableToStructuredGrid1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tableToStructuredGrid1Display_1.OpacityByArray = 0
tableToStructuredGrid1Display_1.OpacityArray = ['POINTS', 'Field 3']
tableToStructuredGrid1Display_1.OpacityArrayComponent = ''
tableToStructuredGrid1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
tableToStructuredGrid1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tableToStructuredGrid1Display_1.SelectionCellLabelBold = 0
tableToStructuredGrid1Display_1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
tableToStructuredGrid1Display_1.SelectionCellLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.SelectionCellLabelFontFile = ''
tableToStructuredGrid1Display_1.SelectionCellLabelFontSize = 18
tableToStructuredGrid1Display_1.SelectionCellLabelItalic = 0
tableToStructuredGrid1Display_1.SelectionCellLabelJustification = 'Left'
tableToStructuredGrid1Display_1.SelectionCellLabelOpacity = 1.0
tableToStructuredGrid1Display_1.SelectionCellLabelShadow = 0
tableToStructuredGrid1Display_1.SelectionPointLabelBold = 0
tableToStructuredGrid1Display_1.SelectionPointLabelColor = [1.0, 1.0, 0.0]
tableToStructuredGrid1Display_1.SelectionPointLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.SelectionPointLabelFontFile = ''
tableToStructuredGrid1Display_1.SelectionPointLabelFontSize = 18
tableToStructuredGrid1Display_1.SelectionPointLabelItalic = 0
tableToStructuredGrid1Display_1.SelectionPointLabelJustification = 'Left'
tableToStructuredGrid1Display_1.SelectionPointLabelOpacity = 1.0
tableToStructuredGrid1Display_1.SelectionPointLabelShadow = 0
tableToStructuredGrid1Display_1.PolarAxes = 'PolarAxesRepresentation'
tableToStructuredGrid1Display_1.ScalarOpacityFunction = None
tableToStructuredGrid1Display_1.ScalarOpacityUnitDistance = 0.06792325611820516
tableToStructuredGrid1Display_1.SelectMapper = 'Projected tetra'
tableToStructuredGrid1Display_1.SamplingDimensions = [128, 128, 128]

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
tableToStructuredGrid1Display_1.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
tableToStructuredGrid1Display_1.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
tableToStructuredGrid1Display_1.GlyphType.TipResolution = 6
tableToStructuredGrid1Display_1.GlyphType.TipRadius = 0.1
tableToStructuredGrid1Display_1.GlyphType.TipLength = 0.35
tableToStructuredGrid1Display_1.GlyphType.ShaftResolution = 6
tableToStructuredGrid1Display_1.GlyphType.ShaftRadius = 0.03
tableToStructuredGrid1Display_1.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToStructuredGrid1Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
tableToStructuredGrid1Display_1.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToStructuredGrid1Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
tableToStructuredGrid1Display_1.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tableToStructuredGrid1Display_1.DataAxesGrid.XTitle = 'X Axis'
tableToStructuredGrid1Display_1.DataAxesGrid.YTitle = 'Y Axis'
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitle = 'Z Axis'
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XTitleOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YTitleOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZTitleOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.FacesToRender = 63
tableToStructuredGrid1Display_1.DataAxesGrid.CullBackface = 0
tableToStructuredGrid1Display_1.DataAxesGrid.CullFrontface = 1
tableToStructuredGrid1Display_1.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.ShowGrid = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ShowEdges = 1
tableToStructuredGrid1Display_1.DataAxesGrid.ShowTicks = 1
tableToStructuredGrid1Display_1.DataAxesGrid.LabelUniqueEdgesOnly = 1
tableToStructuredGrid1Display_1.DataAxesGrid.AxesToLabel = 63
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XLabelOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YLabelOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelFontFile = ''
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelBold = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelItalic = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelFontSize = 12
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelShadow = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZLabelOpacity = 1.0
tableToStructuredGrid1Display_1.DataAxesGrid.XAxisNotation = 'Mixed'
tableToStructuredGrid1Display_1.DataAxesGrid.XAxisPrecision = 2
tableToStructuredGrid1Display_1.DataAxesGrid.XAxisUseCustomLabels = 0
tableToStructuredGrid1Display_1.DataAxesGrid.XAxisLabels = []
tableToStructuredGrid1Display_1.DataAxesGrid.YAxisNotation = 'Mixed'
tableToStructuredGrid1Display_1.DataAxesGrid.YAxisPrecision = 2
tableToStructuredGrid1Display_1.DataAxesGrid.YAxisUseCustomLabels = 0
tableToStructuredGrid1Display_1.DataAxesGrid.YAxisLabels = []
tableToStructuredGrid1Display_1.DataAxesGrid.ZAxisNotation = 'Mixed'
tableToStructuredGrid1Display_1.DataAxesGrid.ZAxisPrecision = 2
tableToStructuredGrid1Display_1.DataAxesGrid.ZAxisUseCustomLabels = 0
tableToStructuredGrid1Display_1.DataAxesGrid.ZAxisLabels = []

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tableToStructuredGrid1Display_1.PolarAxes.Visibility = 0
tableToStructuredGrid1Display_1.PolarAxes.Translation = [0.0, 0.0, 0.0]
tableToStructuredGrid1Display_1.PolarAxes.Scale = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.Orientation = [0.0, 0.0, 0.0]
tableToStructuredGrid1Display_1.PolarAxes.EnableCustomBounds = [0, 0, 0]
tableToStructuredGrid1Display_1.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.EnableCustomRange = 0
tableToStructuredGrid1Display_1.PolarAxes.CustomRange = [0.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.RadialAxesVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.DrawRadialGridlines = 1
tableToStructuredGrid1Display_1.PolarAxes.PolarArcsVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.DrawPolarArcsGridlines = 1
tableToStructuredGrid1Display_1.PolarAxes.NumberOfRadialAxes = 0
tableToStructuredGrid1Display_1.PolarAxes.AutoSubdividePolarAxis = 1
tableToStructuredGrid1Display_1.PolarAxes.NumberOfPolarAxis = 0
tableToStructuredGrid1Display_1.PolarAxes.MinimumRadius = 0.0
tableToStructuredGrid1Display_1.PolarAxes.MinimumAngle = 0.0
tableToStructuredGrid1Display_1.PolarAxes.MaximumAngle = 90.0
tableToStructuredGrid1Display_1.PolarAxes.RadialAxesOriginToPolarAxis = 1
tableToStructuredGrid1Display_1.PolarAxes.Ratio = 1.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitle = 'Radial Distance'
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleLocation = 'Bottom'
tableToStructuredGrid1Display_1.PolarAxes.PolarLabelVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.PolarLabelFormat = '%-#6.3g'
tableToStructuredGrid1Display_1.PolarAxes.PolarLabelExponentLocation = 'Labels'
tableToStructuredGrid1Display_1.PolarAxes.RadialLabelVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.RadialLabelFormat = '%-#3.1f'
tableToStructuredGrid1Display_1.PolarAxes.RadialLabelLocation = 'Bottom'
tableToStructuredGrid1Display_1.PolarAxes.RadialUnitsVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.ScreenSize = 10.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleOpacity = 1.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleFontFile = ''
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleBold = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleItalic = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleShadow = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTitleFontSize = 12
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelOpacity = 1.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelFontFile = ''
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelBold = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelItalic = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelShadow = 0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisLabelFontSize = 12
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextOpacity = 1.0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextFontFile = ''
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextBold = 0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextItalic = 0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextShadow = 0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTextFontSize = 12
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextBold = 0
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextItalic = 0
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextShadow = 0
tableToStructuredGrid1Display_1.PolarAxes.SecondaryRadialAxesTextFontSize = 12
tableToStructuredGrid1Display_1.PolarAxes.EnableDistanceLOD = 1
tableToStructuredGrid1Display_1.PolarAxes.DistanceLODThreshold = 0.7
tableToStructuredGrid1Display_1.PolarAxes.EnableViewAngleLOD = 1
tableToStructuredGrid1Display_1.PolarAxes.ViewAngleLODThreshold = 0.7
tableToStructuredGrid1Display_1.PolarAxes.SmallestVisiblePolarAngle = 0.5
tableToStructuredGrid1Display_1.PolarAxes.PolarTicksVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.ArcTicksOriginToPolarAxis = 1
tableToStructuredGrid1Display_1.PolarAxes.TickLocation = 'Both'
tableToStructuredGrid1Display_1.PolarAxes.AxisTickVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.AxisMinorTickVisibility = 0
tableToStructuredGrid1Display_1.PolarAxes.ArcTickVisibility = 1
tableToStructuredGrid1Display_1.PolarAxes.ArcMinorTickVisibility = 0
tableToStructuredGrid1Display_1.PolarAxes.DeltaAngleMajor = 10.0
tableToStructuredGrid1Display_1.PolarAxes.DeltaAngleMinor = 5.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisMajorTickSize = 0.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTickRatioSize = 0.3
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisMajorTickThickness = 1.0
tableToStructuredGrid1Display_1.PolarAxes.PolarAxisTickRatioThickness = 0.5
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisMajorTickSize = 0.0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTickRatioSize = 0.3
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
tableToStructuredGrid1Display_1.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
tableToStructuredGrid1Display_1.PolarAxes.ArcMajorTickSize = 0.0
tableToStructuredGrid1Display_1.PolarAxes.ArcTickRatioSize = 0.3
tableToStructuredGrid1Display_1.PolarAxes.ArcMajorTickThickness = 1.0
tableToStructuredGrid1Display_1.PolarAxes.ArcTickRatioThickness = 0.5
tableToStructuredGrid1Display_1.PolarAxes.Use2DMode = 0
tableToStructuredGrid1Display_1.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(tableToStructuredGrid1Display_1, ('POINTS', 'Field 3'))

# rescale color and/or opacity maps used to include current data range
tableToStructuredGrid1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tableToStructuredGrid1Display_1.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Field3'
field3LUT = GetColorTransferFunction('Field3')
field3LUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
field3LUT.InterpretValuesAsCategories = 0
field3LUT.AnnotationsInitialized = 0
field3LUT.ShowCategoricalColorsinDataRangeOnly = 0
field3LUT.RescaleOnVisibilityChange = 0
field3LUT.EnableOpacityMapping = 0
field3LUT.RGBPoints = [-1.00117, 0.231373, 0.298039, 0.752941, 0.0, 0.865003, 0.865003, 0.865003, 1.00117, 0.705882, 0.0156863, 0.14902]
field3LUT.UseLogScale = 0
field3LUT.ColorSpace = 'Diverging'
field3LUT.UseBelowRangeColor = 0
field3LUT.BelowRangeColor = [0.0, 0.0, 0.0]
field3LUT.UseAboveRangeColor = 0
field3LUT.AboveRangeColor = [0.5, 0.5, 0.5]
field3LUT.NanColor = [1.0, 1.0, 0.0]
field3LUT.NanOpacity = 1.0
field3LUT.Discretize = 1
field3LUT.NumberOfTableValues = 256
field3LUT.ScalarRangeInitialized = 1.0
field3LUT.HSVWrap = 0
field3LUT.VectorComponent = 0
field3LUT.VectorMode = 'Magnitude'
field3LUT.AllowDuplicateScalars = 1
field3LUT.Annotations = []
field3LUT.ActiveAnnotatedValues = []
field3LUT.IndexedColors = []
field3LUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Field3'
field3PWF = GetOpacityTransferFunction('Field3')
field3PWF.Points = [-1.00117, 0.0, 0.5, 0.0, 1.00117, 1.0, 0.5, 0.0]
field3PWF.AllowDuplicateScalars = 1
field3PWF.UseLogScale = 0
field3PWF.ScalarRangeInitialized = 1

# change representation type
tableToStructuredGrid1Display_1.SetRepresentationType('Surface')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
field3LUT.ApplyPreset('Cold and Hot', True)

# reset view to fit data
renderView1.ResetCamera()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [6.429665074441689, 0.0, 0.0]
renderView1.CameraFocalPoint = [-1e-20, 0.0, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 1.664119774896026

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
