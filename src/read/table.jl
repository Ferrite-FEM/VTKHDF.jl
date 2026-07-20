# Table reader: RowData columns + NumberOfRows (one entry per step, or a
# single entry for static files).

struct ReadTable <: ReaderKind
    rows::Vector{Int}  # per-step row counts (single entry for static files)
end

reader_type_string(::ReadTable) = "Table"
data_locations(::ReadTable) = (VTKRowData(), VTKFieldData())

function read_table_kind(root::HDF5.Group)
    rows = Int.(read(require_dataset(root, "NumberOfRows")))
    all(>=(0), rows) || error("NumberOfRows contains negative counts")
    return ReadTable(rows)
end

function validate_reader(kind::ReadTable, root::HDF5.Group, steps::StepsInfo)
    length(kind.rows) == steps.nsteps ||
        error("NumberOfRows has $(length(kind.rows)) entries, expected NSteps = $(steps.nsteps)")
    return nothing
end

"""
    nrows(r_or_step) -> Int

Number of rows of a Table (for a temporal table: of one time step).
"""
function nrows(r::VTKHDFReader{ReadTable})
    check_open(r)
    is_temporal(r) &&
        error("temporal file: read rows through read_timestep, e.g. nrows(read_timestep(r, i))")
    return isempty(r.kind.rows) ? 0 : r.kind.rows[1]
end
nrows(s::VTKHDFTimeStep{<:VTKHDFReader{ReadTable}}) =
    (check_open(s.reader); s.reader.kind.rows[s.index])

function read_step_array(
        r::VTKHDFReader{ReadTable}, kind::ReadTable,
        loc::VTKRowData, name::String, ds::HDF5.Dataset, i::Int
    )
    n = kind.rows[i]
    default = sum(view(kind.rows, 1:(i - 1)); init = 0)
    off = steps_entry(steps_info(r), "RowDataOffsets/" * name, i, default)
    return read_tuple_rows(ds, (off + 1):(off + n))
end
