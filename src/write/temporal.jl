# Shared temporal ("Steps" group) machinery.
#
# A temporal file appends per-step data to the same datasets used by static
# files and records per-step read offsets in the Steps group. Static geometry
# is stored exactly once: when a step supplies no new geometry, the previous
# step's offsets are repeated (the spec's reuse mechanism).

function steps_group(vtk::VTKHDFFile)
    haskey(vtk.root, "Steps") && return vtk.root["Steps"]::HDF5.Group
    sg = HDF5.create_group(vtk.root, "Steps")
    HDF5.attrs(sg)["NSteps"] = Int64(0)
    return sg
end

"""
    write_timestep(f::Function, vtk, t::Real; geometry...)

Write one time step with time value `t` to a temporal VTKHDF file. Data arrays
are written inside the do-block with the same `vtk[name] = data` syntax used
for static files.

With no geometry keyword arguments, the previous step's geometry is reused and
stored only once. Passing new geometry (`points`/`cells` for unstructured
grids, `x`/`y`/`z` for rectilinear grids, `points` for structured grids)
appends it; this step and later ones then use it.

The set of data arrays (names, element types, component counts) must be the
same for every step; it is fixed by the first step.

```julia
vtk = vtkhdf_grid("sim", points, cells; temporal = true)
for (t, u) in timesteps
    write_timestep(vtk, t) do frame
        frame["u"] = u
    end
end
close(vtk)
```
"""
function write_timestep(f::Function, vtk::VTKHDFFile, t::Real; kwargs...)
    vtk.isopen || error("file is closed")
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    vtk.temporal || error("write_timestep requires a file opened with temporal = true")
    vtk.in_step && error("nested write_timestep")
    vtk.in_step = true
    copy!(vtk.step_start_rows, vtk.data_rows)
    kind_begin_step!(vtk, vtk.kind)
    ok = false
    try
        isempty(kwargs) || supply_step_geometry!(vtk, vtk.kind; kwargs...)
        f(vtk)
        ok = true
    finally
        # A failure anywhere leaves already-appended data unreferenced by any
        # step; mark the file as failed instead of trying to roll back HDF5
        # extents (the data written by completed steps remains valid).
        try
            ok && end_step!(vtk, t)
        catch
            ok = false
            rethrow()
        finally
            vtk.in_step = false
            ok || (vtk.failed = true)
        end
    end
    return vtk
end

kind_begin_step!(vtk, kind::DatasetKind) = nothing

supply_step_geometry!(vtk, kind::DatasetKind; kwargs...) =
    throw(ArgumentError("unsupported geometry keyword arguments $(keys(kwargs)) for $(nameof(typeof(kind)))"))

function end_step!(vtk::VTKHDFFile, t::Real)
    kind = vtk.kind
    geom = finish_step_geometry!(vtk, kind)
    # Fixed array schema across steps, checked before recording the step.
    keys_now = Set(keys(vtk.data_rows))
    if vtk.schema === nothing
        vtk.schema = keys_now
    elseif keys_now != vtk.schema
        extra = sort!(collect(setdiff(keys_now, vtk.schema)))
        missing_ = sort!(collect(setdiff(vtk.schema, keys_now)))
        error(
            "temporal array schema mismatch: new arrays $extra, missing arrays $missing_; " *
                "the array schema is fixed by the first time step"
        )
    end
    validate_data_totals(vtk, step_expected_totals(vtk, kind, geom))
    sg = steps_group(vtk)
    append_rows(appendable(vtk, sg, "Values", Float64, ()), Float64(t))
    append_step_offsets!(vtk, kind, sg, geom)
    append_data_offsets!(vtk, sg)
    vtk.nsteps += 1
    HDF5.attrs(sg)["NSteps"] = Int64(vtk.nsteps)
    push!(vtk.step_values, Float64(t))
    return nothing
end

finish_step_geometry!(vtk, kind::DatasetKind) = nothing
append_step_offsets!(vtk, kind::DatasetKind, sg, geom) = nothing
step_expected_totals(vtk, kind::DatasetKind, geom) = Dict{String, Int}()

# Whether per-array Steps offsets are recorded for point/cell data. Image-like
# kinds use the time-prepended array dimension instead.
uses_data_offsets(kind::DatasetKind) = false

function append_data_offsets!(vtk::VTKHDFFile, sg)
    schema = vtk.schema
    schema === nothing && return nothing
    for key in sort!(collect(schema))
        groupname, name = split(key, '/'; limit = 2)
        start = get(vtk.step_start_rows, key, 0)
        nthis = vtk.data_rows[key] - start
        if groupname == "FieldData"
            nthis > 0 || error("field array $key was not written in this time step")
            og = get_or_create_group(sg, "FieldDataOffsets")
            append_rows(appendable(vtk, og, name, Int64, ()), Int64(start))
            szg = get_or_create_group(sg, "FieldDataSizes")
            ncomp = vtk.field_ncomp[key]
            append_rows(appendable(vtk, szg, name, Int64, (2,)), reshape(Int64[ncomp, nthis], 2, 1))
        elseif uses_data_offsets(vtk.kind)
            offsets_name = groupname * "Offsets"  # PointDataOffsets, CellDataOffsets, RowDataOffsets
            og = get_or_create_group(sg, offsets_name)
            append_rows(appendable(vtk, og, name, Int64, ()), Int64(start))
        end
    end
    return nothing
end
