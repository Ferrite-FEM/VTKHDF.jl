# Read files written by VTK's own vtkHDFWriter (not by this package).
# Gated on the python VTK installation like the writer-validation tests;
# vtkHDFWriter needs VTK >= 9.3.

const HAS_VTK_WRITER = HAS_VTK &&
    success(pipeline(`$PY -c "from vtkmodules.vtkIOHDF import vtkHDFWriter"`; stderr = devnull))

@testset "read files written by VTK" begin
    if !HAS_VTK_WRITER
        @info "Skipping VTK-written interop tests (python VTK with vtkHDFWriter not available)"
    else
        dir = mktempdir()
        run(pipeline(`$PY $(joinpath(@__DIR__, "vtk_write_samples.py")) $dir`; stdout = devnull))

        vtkhdf_open(joinpath(dir, "ug.vtkhdf")) do r
            @test VTKHDF.dataset_type(r) == "UnstructuredGrid"
            @test VTKHDF.npoints(r) == 4
            @test VTKHDF.ncells(r) == 1
            @test read_points(r) == Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
            cells = read_cells(r)
            @test length(cells) == 1
            @test cells[1].ctype == VTKCellTypes.VTK_TETRA
            @test collect(cells[1].connectivity) == [1, 2, 3, 4]
            @test r["u"] == [1.0, 2.0, 3.0, 4.0]
            @test r["c", VTKCellData()] == [7.0]
        end

        vtkhdf_open(joinpath(dir, "pd.vtkhdf")) do r
            @test VTKHDF.dataset_type(r) == "PolyData"
            @test VTKHDF.npoints(r) == 4
            cells = read_cells(r)
            @test length(cells.polygons) == 1
            @test collect(cells.polygons[1].connectivity) == [1, 2, 3, 4]
            @test isempty(cells.lines)
            @test r["u"] == [1.0, 2.0, 3.0, 4.0]
        end

        # written only when this VTK's vtkHDFWriter supports ImageData
        if isfile(joinpath(dir, "img.vtkhdf"))
            vtkhdf_open(joinpath(dir, "img.vtkhdf")) do r
                @test VTKHDF.dataset_type(r) == "ImageData"
                @test VTKHDF.grid_info(r).dims == (2, 3, 1)
                @test VTKHDF.grid_info(r).spacing == (0.5, 1.0, 1.0)
                @test vec(r["u"]) == collect(0.0:5.0)
            end
        else
            @info "vtkHDFWriter does not support ImageData in this VTK; skipping"
        end
    end
end
