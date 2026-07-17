"""Dump the content of a VTKHDF file as JSON using VTK's vtkHDFReader.

Usage: python3 vtk_validate.py FILE

Output: {"time_steps": [...] | null, "steps": [DESCRIPTION, ...]} with one
description per time step (a single one for static files).
"""
import sys
import json

from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonExecutionModel import (
    vtkStreamingDemandDrivenPipeline as sddp,
)


def arrays(dsa):
    out = {}
    for i in range(dsa.GetNumberOfArrays()):
        name = dsa.GetArrayName(i)
        arr = dsa.GetArray(i) if hasattr(dsa, "GetArray") else None
        if arr is None:
            a = dsa.GetAbstractArray(i)
            out[name] = [a.GetValue(j) for j in range(a.GetNumberOfValues())]
        else:
            out[name] = vtk_to_numpy(arr).tolist()
    return out


def describe_leaf(ds):
    d = {"class": ds.GetClassName()}
    for attr, key in (
        ("GetNumberOfPoints", "npoints"),
        ("GetNumberOfCells", "ncells"),
        ("GetNumberOfRows", "nrows"),
    ):
        if hasattr(ds, attr):
            d[key] = getattr(ds, attr)()
    if hasattr(ds, "GetPoints") and callable(getattr(ds, "GetPoints", None)):
        try:
            pts = ds.GetPoints()
            if pts is not None and hasattr(pts, "GetData"):
                d["points"] = vtk_to_numpy(pts.GetData()).tolist()
        except TypeError:  # e.g. vtkImageData.GetPoints needs arguments
            pass
    if hasattr(ds, "GetDimensions"):
        try:
            d["dimensions"] = list(ds.GetDimensions())
        except TypeError:
            pass
    if hasattr(ds, "GetOrigin"):
        d["origin"] = list(ds.GetOrigin())
    if hasattr(ds, "GetSpacing"):
        try:
            d["spacing"] = list(ds.GetSpacing())
        except TypeError:
            pass
    if ds.GetClassName() == "vtkImageData":
        m = ds.GetDirectionMatrix()
        d["direction"] = [m.GetElement(i, j) for i in range(3) for j in range(3)]
    if hasattr(ds, "GetCellType"):
        d["cell_types"] = [ds.GetCellType(i) for i in range(ds.GetNumberOfCells())]
        conn = []
        for i in range(ds.GetNumberOfCells()):
            c = ds.GetCell(i)
            conn.append([c.GetPointId(j) for j in range(c.GetNumberOfPoints())])
        d["connectivity"] = conn
    if hasattr(ds, "GetXCoordinates"):
        for f, key in (
            ("GetXCoordinates", "x"),
            ("GetYCoordinates", "y"),
            ("GetZCoordinates", "z"),
        ):
            d[key] = vtk_to_numpy(getattr(ds, f)()).tolist()
    for get, key in (
        ("GetPointData", "point_data"),
        ("GetCellData", "cell_data"),
        ("GetFieldData", "field_data"),
        ("GetRowData", "row_data"),
    ):
        if hasattr(ds, get):
            d[key] = arrays(getattr(ds, get)())
    return d


def describe(obj):
    cls = obj.GetClassName()
    if cls == "vtkPartitionedDataSetCollection":
        d = {"class": cls, "blocks": []}
        for i in range(obj.GetNumberOfPartitionedDataSets()):
            pds = obj.GetPartitionedDataSet(i)
            parts = [describe(pds.GetPartitionAsDataObject(j))
                     for j in range(pds.GetNumberOfPartitions())]
            d["blocks"].append(parts[0] if len(parts) == 1 else {"class": "vtkPartitionedDataSet", "partitions": parts})
        asm = obj.GetDataAssembly()
        if asm is not None:
            def node(n):
                return {
                    "name": asm.GetNodeName(n),
                    "children": [node(c) for c in asm.GetChildNodes(n, False)],
                    "datasets": list(asm.GetDataSetIndices(n, False)),
                }
            d["assembly"] = node(0)
        return d
    if cls == "vtkMultiBlockDataSet":
        return {"class": cls,
                "blocks": [describe(obj.GetBlock(i)) for i in range(obj.GetNumberOfBlocks())]}
    if cls == "vtkPartitionedDataSet":
        return {"class": cls,
                "partitions": [describe(obj.GetPartitionAsDataObject(i))
                               for i in range(obj.GetNumberOfPartitions())]}
    if cls == "vtkOverlappingAMR":
        levels = []
        for lvl in range(obj.GetNumberOfLevels()):
            try:
                n = obj.GetNumberOfBlocks(lvl)
            except (AttributeError, TypeError):
                n = obj.GetNumberOfDataSets(lvl)
            spacing = [0.0, 0.0, 0.0]
            obj.GetSpacing(lvl, spacing)
            levels.append({
                "spacing": spacing,
                "boxes": [describe_leaf(obj.GetDataSet(lvl, i)) for i in range(n)],
            })
        bounds = [0.0] * 6
        obj.GetBounds(bounds)
        return {"class": cls, "bounds": bounds, "levels": levels}
    return describe_leaf(obj)


def main():
    fn = sys.argv[1]
    reader = vtkHDFReader()
    reader.SetFileName(fn)
    reader.UpdateInformation()
    info = reader.GetOutputInformation(0)
    ts = info.Get(sddp.TIME_STEPS()) if info.Has(sddp.TIME_STEPS()) else None
    result = {"time_steps": list(ts) if ts else None, "steps": []}
    if ts:
        # A fresh reader per step: VTK 9.6's vtkHDFReader caches composite
        # temporal data and returns stale steps otherwise (also for files
        # written by vtkHDFWriter itself).
        for t in ts:
            reader = vtkHDFReader()
            reader.SetFileName(fn)
            reader.UpdateTimeStep(t)
            result["steps"].append(describe(reader.GetOutput()))
    else:
        reader.Update()
        result["steps"].append(describe(reader.GetOutput()))
    json.dump(result, sys.stdout)


if __name__ == "__main__":
    main()
