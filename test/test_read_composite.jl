# Reading back composite files (PartitionedDataSetCollection /
# MultiBlockDataSet), including the Assembly tree and block lifetimes.

@testset "read collection" begin
    dir = mktempdir()
    fn = joinpath(dir, "col")
    pts = Float64[0 1 0; 0 0 1; 0 0 0]
    cls = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 2, 3])]
    Timg = [1.0 2.0; 3.0 4.0]
    vtkhdf_collection(fn) do col
        mesh = vtkhdf_grid(col, "mesh", pts, cls; temporal = true)
        write_timestep(f -> f["u"] = [1.0, 2.0, 3.0], mesh, 0.0)
        write_timestep(f -> f["u"] = [4.0, 5.0, 6.0], mesh, 1.0)
        img = vtkhdf_grid(col, "img", VTKImageData(), (2, 2))
        img["T"] = Timg
        e = add_empty_block(col, "placeholder")
        n1 = add_node(col, "group1")
        add_block_ref(n1, mesh)
        n2 = add_node(n1, "inner")
        add_block_ref(n2, img)
        add_block_ref(col, e)
        add_block_ref(n2, mesh)  # a block may be referenced repeatedly
    end

    col = vtkhdf_open(fn)
    @test col isa VTKHDF.VTKHDFCollectionReader
    @test VTKHDF.dataset_type(col) == "PartitionedDataSetCollection"
    @test VTKHDF.file_version(col) == (2, 1)
    @test keys(col) == ["mesh", "img", "placeholder"]
    @test haskey(col, "mesh")
    @test !haskey(col, "nope")
    @test_throws ErrorException col["nope"]
    @test occursin("3 blocks", sprint(show, col))

    mesh = col["mesh"]
    @test mesh === col["mesh"]  # block readers are cached
    @test VTKHDF.dataset_type(mesh) == "UnstructuredGrid"
    @test VTKHDF.is_temporal(mesh)
    @test read_points(read_timestep(mesh, 1)) == pts
    @test read_timestep(mesh, 2)["u"] == [4.0, 5.0, 6.0]

    img = col["img"]
    @test img["T"] == reshape(Timg, 2, 2, 1)
    @test VTKHDF.dataset_type(col["placeholder"]) === nothing

    @test VTKHDF.is_temporal(col)
    @test VTKHDF.nsteps(col) == 2
    @test VTKHDF.time_values(col) == [0.0, 1.0]

    asm = VTKHDF.read_assembly(col)
    @test asm.name == "Assembly"
    @test asm.blocks == ["placeholder"]
    @test length(asm.children) == 1
    g1 = asm.children[1]
    @test g1.name == "group1"
    @test g1.blocks == ["mesh"]
    g2 = only(g1.children)
    @test g2.name == "inner"
    @test g2.blocks == ["img", "mesh"]

    # closing the collection invalidates all block readers
    empty_block = col["placeholder"]
    close(col)
    @test !isopen(col)
    @test_throws ErrorException keys(col)
    @test_throws ErrorException VTKHDF.dataset_type(col)
    @test_throws ErrorException read_timestep(mesh, 1)
    @test_throws ErrorException VTKHDF.dataset_type(mesh)
    @test_throws ErrorException VTKHDF.dataset_type(empty_block)
    @test_throws ErrorException img["T"]
    close(col)  # no-op
end

@testset "read multiblock" begin
    dir = mktempdir()
    fn = joinpath(dir, "mb")
    vtkhdf_multiblock(fn) do col
        t = vtkhdf_table(col, "tbl")
        t["a"] = [1.0]
        vtkhdf_grid(col, "grid", VTKImageData(), (2, 2))
    end
    vtkhdf_open(fn) do col
        @test VTKHDF.dataset_type(col) == "MultiBlockDataSet"
        @test keys(col) == ["tbl", "grid"]
        @test col["tbl"]["a"] == [1.0]
        @test !VTKHDF.is_temporal(col)
        @test_throws ErrorException VTKHDF.nsteps(col)
    end

    # do-block form returns the collection
    col = vtkhdf_open(identity, fn)
    @test !isopen(col)
end

@testset "read assembly adversarial links" begin
    dir = mktempdir()
    # A soft link inside the Assembly pointing at another assembly node makes
    # a node resolve twice; the reader must reject it rather than build a
    # bogus tree.
    fn = joinpath(dir, "badasm.vtkhdf")
    vtkhdf_collection(fn) do col
        add_node(col, "a")
        add_node(col, "b")
    end
    h5open(fn, "r+") do file
        VTKHDF.create_soft_link(file["VTKHDF/Assembly/a"]::HDF5.Group, "sneaky", "/VTKHDF/Assembly/b")
    end
    vtkhdf_open(fn) do col
        @test_throws ErrorException VTKHDF.read_assembly(col)
    end
end

@testset "read collection temporal consistency" begin
    dir = mktempdir()
    fn = joinpath(dir, "cols")
    vtkhdf_collection(fn) do col
        a = vtkhdf_table(col, "a"; temporal = true)
        write_timestep(f -> f["x"] = [1.0], a, 0.0)
        b = vtkhdf_table(col, "b"; temporal = true)
        write_timestep(f -> f["y"] = [2.0], b, 0.0)
        vtkhdf_table(col, "static")["z"] = [3.0]
    end
    vtkhdf_open(fn) do col
        @test VTKHDF.is_temporal(col)
        @test VTKHDF.nsteps(col) == 1
        @test !VTKHDF.is_temporal(col["static"])
        @test col["static"]["z"] == [3.0]
        @test read_timestep(col["a"], 1)["x"] == [1.0]
    end
end
