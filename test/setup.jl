# Shared test helpers, loaded into every test sandbox module.
using VTKHDF
using Test
using HDF5
import JSON

const PY = Sys.which("python3")
const HAS_VTK = PY !== nothing &&
    success(pipeline(`$PY -c "import vtkmodules.vtkIOHDF"`; stderr = devnull))
const HAS_H5DUMP = Sys.which("h5dump") !== nothing
const VALIDATE = joinpath(@__DIR__, "vtk_validate.py")

# On CI, at least one job must run the authoritative VTK/h5dump validation;
# it sets this variable so that quietly skipping is impossible.
if get(ENV, "VTKHDF_REQUIRE_VTK", "") == "true"
    HAS_VTK || error("VTKHDF_REQUIRE_VTK is set but python3 + vtkmodules is not available")
    HAS_H5DUMP || error("VTKHDF_REQUIRE_VTK is set but h5dump is not available")
    # the reader-interop tests write reference files with vtkHDFWriter
    # (VTK >= 9.3); they must not silently skip on the validating CI job
    success(pipeline(`$PY -c "from vtkmodules.vtkIOHDF import vtkHDFWriter"`; stderr = devnull)) ||
        error("VTKHDF_REQUIRE_VTK is set but this VTK has no vtkHDFWriter (VTK >= 9.3 required)")
end

# Read a file back through VTK's own reader as a JSON-compatible structure.
vtkdump(path) = JSON.parse(read(pipeline(`$PY $VALIDATE $path`; stderr = devnull), String))

# Raw on-disk dataspace dims via h5dump (independent of HDF5.jl's dimension
# reversal, so layout bugs cannot cancel out).
function h5dump_shape(path, dataset)
    out = read(`h5dump -d $dataset --header $path`, String)
    m = match(r"DATASPACE\s+SIMPLE\s*\{\s*\(\s*([0-9, ]+)\)", out)
    m === nothing && (m = match(r"DATASPACE\s+SCALAR", out); return ())
    return Tuple(parse.(Int, split(m.captures[1], ',')))
end

nested(points::AbstractMatrix) = [points[:, i] for i in axes(points, 2)]
