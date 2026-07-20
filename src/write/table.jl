# Table writer (spec 2.8): NumberOfRows (one entry per step) + RowData columns.

"""
    VTKRowData()

The data location of Table columns, usable wherever `VTKPointData()` etc.
are: `tbl["a", VTKRowData()] = col` when writing, `r["a", VTKRowData()]` and
`keys(r, VTKRowData())` when reading.
"""
struct VTKRowData <: VTKBase.AbstractFieldData end
location_group(::VTKRowData) = "RowData"

mutable struct TableState <: DatasetKind
    cum_rows::Int
end

"""
    vtkhdf_table(filename; temporal = false, compress = false, kwargs...)

Create a VTKHDF `Table` file. Columns are written with `tbl[name] = column`
where a column is a vector or a `(ncomponents, nrows)` matrix; all columns
must have the same number of rows. For temporal tables, write columns inside
[`write_timestep`](@ref).

Supports the do-block form `vtkhdf_table(fname) do tbl ... end`.
"""
function vtkhdf_table(filename::AbstractString; kwargs...)
    return init_table(filename; kwargs...)
end

function init_table(dest; temporal = false, compress = false, chunk_size = 0)
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "Table")
    write_version_attribute(root, (2, 8))
    return make_vtkfile(
        file, root, TableState(0);
        temporal, compress, chunk_size, version = (2, 8)
    )
end

resolve_location(vtk, kind::TableState, data) = VTKRowData()

write_array!(vtk::VTKHDFFile, kind::TableState, loc::VTKRowData, name::AbstractString, data) =
    append_tuple_data!(vtk, "RowData", name, data)

write_array!(vtk::VTKHDFFile, kind::TableState, loc::Union{VTKPointData, VTKCellData}, name::AbstractString, data) =
    throw(ArgumentError("Table files only have row data; write columns with tbl[name] = data"))

# Consistent row count over the given group entries (delta between two
# data_rows snapshots).
function consistent_rows(rows::Dict{String, Int}, start::Dict{String, Int}, groupname::String)
    n = -1
    for (key, total) in rows
        startswith(key, groupname * "/") || continue
        delta = total - get(start, key, 0)
        if n == -1
            n = delta
        elseif n != delta
            error("table columns have inconsistent row counts ($n vs $delta)")
        end
    end
    return max(n, 0)
end

function finalize_kind!(vtk, kind::TableState)
    get_or_create_group(vtk.root, "RowData")
    if vtk.temporal
        # zero-step file: materialize an empty NumberOfRows
        vtk.nsteps == 0 && appendable(vtk, vtk.root, "NumberOfRows", Int64, ())
        return nothing
    end
    nrows = consistent_rows(vtk.data_rows, Dict{String, Int}(), "RowData")
    HDF5.write_dataset(vtk.root, "NumberOfRows", Int64[nrows])
    return nothing
end

# ---- temporal hooks ----

uses_data_offsets(::TableState) = true

function finish_step_geometry!(vtk, kind::TableState)
    nrows = consistent_rows(vtk.data_rows, vtk.step_start_rows, "RowData")
    kind.cum_rows += nrows
    return nrows
end

step_expected_totals(vtk, kind::TableState, nrows) = Dict("RowData" => kind.cum_rows)

function append_step_offsets!(vtk, kind::TableState, sg, nrows)
    append_rows(appendable(vtk, vtk.root, "NumberOfRows", Int64, ()), Int64(nrows))
    return nothing
end
