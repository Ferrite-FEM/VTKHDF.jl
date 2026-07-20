"""Write small reference VTKHDF files with VTK's own vtkHDFWriter.

Usage: python3 vtk_write_samples.py OUTDIR

Writes ug.vtkhdf (UnstructuredGrid) and pd.vtkhdf (PolyData) with known
contents; the Julia test reads them back and compares.
"""
import sys

from vtkmodules.vtkCommonCore import vtkPoints, vtkDoubleArray
from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
    vtkPolyData,
    vtkCellArray,
    VTK_TETRA,
)
from vtkmodules.vtkIOHDF import vtkHDFWriter

outdir = sys.argv[1]


def named_array(name, values):
    arr = vtkDoubleArray()
    arr.SetName(name)
    for v in values:
        arr.InsertNextValue(v)
    return arr


# UnstructuredGrid: one tetrahedron
ug = vtkUnstructuredGrid()
pts = vtkPoints()
for p in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]:
    pts.InsertNextPoint(p)
ug.SetPoints(pts)
ug.Allocate(1)
ug.InsertNextCell(VTK_TETRA, 4, [0, 1, 2, 3])
ug.GetPointData().AddArray(named_array("u", [1.0, 2.0, 3.0, 4.0]))
ug.GetCellData().AddArray(named_array("c", [7.0]))

w = vtkHDFWriter()
w.SetFileName(outdir + "/ug.vtkhdf")
w.SetInputData(ug)
w.Write()

# PolyData: a quad polygon
pd = vtkPolyData()
pts = vtkPoints()
for p in [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]:
    pts.InsertNextPoint(p)
pd.SetPoints(pts)
polys = vtkCellArray()
polys.InsertNextCell(4, [0, 1, 2, 3])
pd.SetPolys(polys)
pd.GetPointData().AddArray(named_array("u", [1.0, 2.0, 3.0, 4.0]))

w = vtkHDFWriter()
w.SetFileName(outdir + "/pd.vtkhdf")
w.SetInputData(pd)
w.Write()

# ImageData support in vtkHDFWriter is newer than UnstructuredGrid/PolyData
# (not present up to at least VTK 9.6); write it best-effort — with VTK's
# error chatter silenced — and let the Julia test skip when absent.
try:
    from vtkmodules.vtkCommonCore import vtkObject
    from vtkmodules.vtkCommonDataModel import vtkImageData

    vtkObject.GlobalWarningDisplayOff()
    img = vtkImageData()
    img.SetDimensions(2, 3, 1)
    img.SetSpacing(0.5, 1.0, 1.0)
    img.GetPointData().AddArray(named_array("u", [float(i) for i in range(6)]))
    w = vtkHDFWriter()
    w.SetFileName(outdir + "/img.vtkhdf")
    w.SetInputData(img)
    if w.Write() != 1:
        import os

        if os.path.exists(outdir + "/img.vtkhdf"):
            os.remove(outdir + "/img.vtkhdf")
except Exception:
    pass

print("ok")
