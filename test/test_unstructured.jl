@testset "UnstructuredGrid" begin
    dir = mktempdir()
    points = Float64[
        0 1 1 0 0 1 1 0 0.5
        0 0 1 1 0 0 1 1 0.5
        0 0 0 0 1 1 1 1 2.0
    ]
    cells = [
        MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8),
        MeshCell(VTKCellTypes.VTK_TETRA, [5, 6, 7, 9]),
    ]
    temp = Float64.(1:9)
    vel = rand(3, 9)
    mat = Int32[1, 2]

    @testset "static" begin
        fn = joinpath(dir, "ug.vtkhdf")
        vtkhdf_grid(fn, points, cells) do vtk
            vtk["temp"] = temp
            vtk["vel", VTKPointData(), attribute = :Vectors] = vel
            vtk["mat"] = mat
            vtk["title", VTKFieldData()] = "static file"
            vtk["ke", VTKFieldData()] = [1.5, 2.5]
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "UnstructuredGrid"
            @test attrs(g)["Version"] == [2, 6]  # bumped by attribute marking
            @test read(g["NumberOfPoints"]) == [9]
            @test read(g["NumberOfCells"]) == [2]
            @test read(g["NumberOfConnectivityIds"]) == [12]
            @test read(g["Connectivity"]) == [0:7; [4, 5, 6, 8]]
            @test read(g["Offsets"]) == [0, 8, 12]
            @test read(g["Types"]) == UInt8[12, 10]
            @test read(g["Points"]) == points
            @test read(g["PointData/temp"]) == temp
            @test read(g["PointData/vel"]) == vel
            @test eltype(read(g["CellData/mat"])) == Int32
            @test attrs(g["PointData"])["Vectors"] == "vel"
            @test attrs(g["PointData/vel"])["Attribute"] == "Vectors"
            @test read(g["FieldData/title"]) == ["static file"]
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn, "/VTKHDF/Points") == (9, 3)
            @test h5dump_shape(fn, "/VTKHDF/PointData/vel") == (9, 3)
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test d["time_steps"] === nothing
            s = only(d["steps"])
            @test s["class"] == "vtkUnstructuredGrid"
            @test s["npoints"] == 9 && s["ncells"] == 2
            @test s["cell_types"] == [12, 10]
            @test s["connectivity"] == [collect(0:7), [4, 5, 6, 8]]
            @test s["points"] == nested(points)
            @test s["point_data"]["temp"] == temp
            @test s["point_data"]["vel"] == nested(vel)
            @test s["cell_data"]["mat"] == mat
            @test s["field_data"]["ke"] == [1.5, 2.5]
            @test s["field_data"]["title"] == ["static file"]
        end
    end

    @testset "point input forms" begin
        for (name, pts) in (
                ("tuple", (points[1, :], points[2, :], points[3, :])),
                ("vectors", [Tuple(points[:, i]) for i in axes(points, 2)]),
            )
            fn = joinpath(dir, "ug_$name.vtkhdf")
            vtkhdf_grid(identity, fn, pts, cells)
            h5open(fn) do f
                @test read(f["VTKHDF/Points"]) == points
            end
        end
        # 2-D points are zero-padded
        fn = joinpath(dir, "ug_2d.vtkhdf")
        tri = [MeshCell(VTKCellTypes.VTK_TRIANGLE, 1:3)]
        vtkhdf_grid(identity, fn, Float64[0 1 0; 0 0 1], tri)
        h5open(fn) do f
            @test read(f["VTKHDF/Points"])[3, :] == zeros(3)
        end
        # abstractly typed cell vector (as produced by e.g. `MeshCell[...]`)
        fn = joinpath(dir, "ug_loose.vtkhdf")
        vtkhdf_grid(identity, fn, points, MeshCell[cells...])
        h5open(fn) do f
            @test read(f["VTKHDF/Types"]) == UInt8[12, 10]
        end
    end

    @testset "temporal static geometry" begin
        fn = joinpath(dir, "ug_t.vtkhdf")
        nsteps = 4
        us = [rand(9) for _ in 1:nsteps]
        vtk = vtkhdf_grid(fn, points, cells; temporal = true)
        for s in 1:nsteps
            write_timestep(vtk, 0.5 * (s - 1)) do frame
                frame["u"] = us[s]
                frame["c", VTKCellData()] = Float64[s, s + 1]
                frame["time", VTKFieldData()] = [0.5 * (s - 1)]
            end
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            # geometry stored exactly once
            @test size(g["Points"]) == (3, 9)
            @test length(g["Connectivity"]) == 12
            @test read(g["NumberOfPoints"]) == [9]
            @test length(g["PointData/u"]) == 9 * nsteps
            @test attrs(g["Steps"])["NSteps"] == nsteps
            @test read(g["Steps/Values"]) == 0.5 .* (0:(nsteps - 1))
            @test read(g["Steps/PointOffsets"]) == zeros(nsteps)
            @test read(g["Steps/PartOffsets"]) == zeros(nsteps)
            @test read(g["Steps/NumberOfParts"]) == ones(nsteps)
            @test read(g["Steps/CellOffsets"]) == zeros(1, nsteps)
            @test read(g["Steps/PointDataOffsets/u"]) == 9 .* (0:(nsteps - 1))
            @test read(g["Steps/CellDataOffsets/c"]) == 2 .* (0:(nsteps - 1))
            @test read(g["Steps/FieldDataOffsets/time"]) == 0:(nsteps - 1)
            @test read(g["Steps/FieldDataSizes/time"]) == repeat([1, 1], 1, nsteps)
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test d["time_steps"] == 0.5 .* (0:(nsteps - 1))
            for s in 1:nsteps
                st = d["steps"][s]
                @test st["npoints"] == 9 && st["ncells"] == 2
                @test st["points"] == nested(points)
                @test st["point_data"]["u"] ≈ us[s]
                @test st["cell_data"]["c"] == [s, s + 1]
            end
        end
    end

    @testset "temporal changing geometry" begin
        fn = joinpath(dir, "ug_tc.vtkhdf")
        vtk = vtkhdf_grid(VTKUnstructuredGrid(), fn; temporal = true)
        hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]
        for s in 1:3
            pts = points[:, 1:8] .+ (s - 1.0)
            write_timestep(vtk, Float64(s); points = pts, cells = hex) do frame
                frame["u"] = fill(Float64(s), 8)
            end
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test size(g["Points"]) == (3, 24)
            @test read(g["NumberOfPoints"]) == [8, 8, 8]
            @test read(g["Steps/PointOffsets"]) == [0, 8, 16]
            @test read(g["Steps/PartOffsets"]) == [0, 1, 2]
            @test read(g["Steps/CellOffsets"]) == reshape([0, 1, 2], 1, 3)
            @test read(g["Steps/ConnectivityIdOffsets"]) == reshape([0, 8, 16], 1, 3)
        end
        if HAS_VTK
            d = vtkdump(fn)
            for s in 1:3
                st = d["steps"][s]
                @test st["points"][1][1] == s - 1.0
                @test st["point_data"]["u"] == fill(s, 8)
            end
        end
    end

    @testset "partitions" begin
        fn = joinpath(dir, "ug_parts.vtkhdf")
        cube = points[:, 1:8]
        hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]
        vtk = vtkhdf_grid(VTKUnstructuredGrid(), fn)
        add_partition(
            vtk, cube, hex;
            pointdata = ("u" => Float64.(1:8),), celldata = ("c" => [1.0],)
        )
        add_partition(
            vtk, cube .+ 3.0, hex;
            pointdata = ("u" => Float64.(9:16),), celldata = ("c" => [2.0],)
        )
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["NumberOfPoints"]) == [8, 8]
            @test read(g["Offsets"]) == [0, 8, 0, 8]
            @test read(g["PointData/u"]) == 1:16
        end
        if HAS_VTK
            d = vtkdump(fn)
            s = only(d["steps"])
            @test s["class"] == "vtkPartitionedDataSet"
            @test [p["point_data"]["u"] for p in s["partitions"]] == [1:8, 9:16]
        end
    end

    @testset "polyhedra" begin
        fn = joinpath(dir, "ug_poly.vtkhdf")
        cube = points[:, 1:8]
        poly = VTKPolyhedron(
            1:8,
            (1, 4, 3, 2), (1, 5, 8, 4), (5, 6, 7, 8),
            (6, 2, 3, 7), (1, 2, 6, 5), (3, 4, 8, 7)
        )
        # regular cell first => polyhedron datasets are backfilled
        vtk = vtkhdf_grid(VTKUnstructuredGrid(), fn)
        add_partition(vtk, cube, [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)])
        add_partition(vtk, cube .+ 3.0, [poly])
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Version"] == [2, 5]
            @test read(g["NumberOfFaces"]) == [0, 6]
            @test read(g["NumberOfPolyhedronToFaceIds"]) == [0, 6]
            @test read(g["NumberOfFaceConnectivityIds"]) == [0, 24]
            @test read(g["FaceOffsets"]) == [0; 0:4:24]
            @test read(g["PolyhedronOffsets"]) == [0, 0, 0, 6]
            @test read(g["PolyhedronToFaces"]) == 0:5
            @test read(g["FaceConnectivity"])[1:4] == [0, 3, 2, 1]
        end
        if HAS_VTK
            d = vtkdump(fn)
            parts = only(d["steps"])["partitions"]
            @test parts[2]["cell_types"] == [42]
        end
    end

    @testset "temporal polyhedra" begin
        fn = joinpath(dir, "ug_tpoly.vtkhdf")
        cube = points[:, 1:8]
        poly = VTKPolyhedron(
            1:8,
            (1, 4, 3, 2), (1, 5, 8, 4), (5, 6, 7, 8),
            (6, 2, 3, 7), (1, 2, 6, 5), (3, 4, 8, 7)
        )
        vtk = vtkhdf_grid(VTKUnstructuredGrid(), fn; temporal = true)
        for s in 1:2
            write_timestep(vtk, Float64(s); points = cube .+ s, cells = [poly]) do frame
                frame["u"] = fill(Float64(s), 8)
            end
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["Steps/FaceConnectivityOffsets"]) == [0, 24]
            @test read(g["Steps/FaceOffsetsOffsets"]) == [0, 6]
            @test read(g["Steps/PolyhedronToFaceIdOffsets"]) == [0, 6]
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test [st["cell_types"] for st in d["steps"]] == [[42], [42]]
            @test d["steps"][2]["point_data"]["u"] == fill(2.0, 8)
        end
    end

    @testset "compression" begin
        fn = joinpath(dir, "ug_comp.vtkhdf")
        vtkhdf_grid(identity, fn, points, cells; compress = true)
        h5open(fn) do f
            filters = HDF5.get_create_properties(f["VTKHDF/Points"]).filters
            @test any(f -> f isa HDF5.Filters.Deflate, filters)
        end
        HAS_VTK && @test only(vtkdump(fn)["steps"])["npoints"] == 9
    end
end
