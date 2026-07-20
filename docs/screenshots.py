# Renders the example screenshots in docs/src/assets/examples with ParaView.
# Each scene is saved twice on a transparent background: <name>-light.png with
# dark annotations and <name>-dark.png with light ones; docs/src/assets/custom.css
# shows the right one for the active Documenter theme.
#
# Usage (from the repo root):
#   julia --project -e 'foreach(ex -> include(joinpath("examples", ex)), readdir("examples"))'
#   pvbatch docs/screenshots.py . docs/src/assets/examples
import sys

from paraview.simple import *

datadir, outdir = sys.argv[1], sys.argv[2]
RES = [1200, 900]
BARS = []  # scalar bars whose text flips between the variants
INK = []  # (display, light_rgb, dark_rgb) solid-colored props that flip


def new_view():
    view = CreateRenderView()
    SetActiveView(view)
    view.ViewSize = RES
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    return view


def ink(display, light, dark):
    INK.append((display, light, dark))


def finish(view, name, azimuth=30, elevation=25, zoom=1.0):
    view.ResetCamera(False)
    cam = GetActiveCamera()
    cam.Azimuth(azimuth)
    cam.Elevation(elevation)
    view.ResetCamera(False)
    cam.Zoom(zoom)
    for variant, text in (("light", [0.0, 0.0, 0.0]), ("dark", [0.9, 0.9, 0.9])):
        for bar in BARS:
            bar.TitleColor = text
            bar.LabelColor = text
        for display, light, dark in INK:
            color = light if variant == "light" else dark
            display.AmbientColor = color
            display.DiffuseColor = color
        Render()
        SaveScreenshot(
            outdir + "/" + name + "-" + variant + ".png", view,
            ImageResolution=RES, TransparentBackground=1,
        )
    BARS.clear()
    INK.clear()
    Delete(view)


def colorbar(display, view, array, title=None):
    ColorBy(display, array)
    display.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(array[1])
    bar = GetScalarBar(lut, view)
    bar.Title = title or array[1]
    bar.ComponentTitle = ""
    bar.WindowLocation = "Lower Right Corner"
    bar.TitleFontSize = 22
    bar.LabelFontSize = 20
    bar.TitleBold = 1
    bar.LabelBold = 1
    BARS.append(bar)
    return lut


# --- ImageData: mandelbrot volume rendering -------------------------------
view = new_view()
r = OpenDataFile(datadir + "/mandelbrot.vtkhdf")
d = Show(r, view)
d.SetRepresentationType("Volume")
lut = colorbar(d, view, ("POINTS", "Iterations"))
lut.ApplyPreset("Inferno", True)
finish(view, "image_data", azimuth=40, elevation=25)

# --- UnstructuredGrid: can partitions -------------------------------------
view = new_view()
r = OpenDataFile(datadir + "/can.vtkhdf")
d = Show(r, view)
d.Representation = "Surface With Edges"
colorbar(d, view, ("CELLS", "EQPS"))
g = Glyph(
    Input=r,
    GlyphType="Arrow",
    OrientationArray=["POINTS", "VEL"],
    ScaleArray=["POINTS", "VEL"],
    ScaleFactor=0.2,
    GlyphMode="Every Nth Point",
    Stride=9,
)
gd = Show(g, view)
ink(gd, [0.35, 0.35, 0.35], [0.7, 0.7, 0.7])
finish(view, "unstructured_grid", azimuth=30, elevation=25)

# --- PolyData: warped torus ------------------------------------------------
view = new_view()
r = OpenDataFile(datadir + "/torus.vtkhdf")
w = WarpByVector(Input=r, Vectors=["POINTS", "Warping"])
d = Show(w, view)
d.Representation = "Surface With Edges"
lut = colorbar(d, view, ("CELLS", "Materials"))
lut.InterpretValuesAsCategories = 1
lut.Annotations = ["1", "Material 1", "2", "Material 2"]
lut.IndexedColors = [0.21, 0.47, 0.53, 0.91, 0.59, 0.31]
finish(view, "poly_data", azimuth=20, elevation=55, zoom=1.15)

# --- OverlappingAMR: gaussian pulse ----------------------------------------
view = new_view()
r = OpenDataFile(datadir + "/gaussian_pulse.vtkhdf")
clip = Clip(Input=r, ClipType="Plane")
clip.ClipType.Origin = [-0.75, -0.75, 1.25]
clip.ClipType.Normal = [0, 1, 0]
clip.Invert = 1
d = Show(clip, view)
d.Representation = "Surface With Edges"
lut = colorbar(d, view, ("CELLS", "Gaussian-Pulse"))
lut.ApplyPreset("Viridis", True)
finish(view, "overlapping_amr", azimuth=30, elevation=30)

# --- PartitionedDataSetCollection: composite -------------------------------
view = new_view()
r = OpenDataFile(datadir + "/composite.vtkhdf")
solid = ExtractBlock(Input=r, Assembly="Hierarchy", Selectors=["/Root/Solid"])
ds = Show(solid, view)
ds.Representation = "Surface With Edges"
lut = colorbar(ds, view, ("POINTS", "Temperature"))
surface = ExtractBlock(Input=r, Assembly="Hierarchy", Selectors=["/Root/Surface"])
t = Transform(Input=surface)
t.Transform.Translate = [1.4, 0, 0]
dt = Show(t, view)
dt.Representation = "Surface With Edges"
ColorBy(dt, ("POINTS", None))
dt.AmbientColor = [0.91, 0.59, 0.31]
dt.DiffuseColor = [0.91, 0.59, 0.31]
finish(view, "partitioned_collection", azimuth=30, elevation=25)

# --- Temporal PolyData: travelling wave ------------------------------------
view = new_view()
r = OpenDataFile(datadir + "/wave.vtkhdf")
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
view.ViewTime = 0.3
w = WarpByVector(Input=r, Vectors=["POINTS", "displacement"])
d = Show(w, view)
d.Representation = "Surface With Edges"
lut = colorbar(d, view, ("POINTS", "height"))
lut.ApplyPreset("Cool to Warm (Extended)", True)
finish(view, "temporal_poly_data", azimuth=40, elevation=30, zoom=1.0)

print("all screenshots written to", outdir)
