# Conversion of VTKBase cells to VTKHDF connectivity arrays.

const AnyCell = Union{MeshCell, VTKPolyhedron}

cell_type_id(cell::MeshCell{VTKCellType}) = cell.ctype.vtk_id % UInt8
cell_type_id(::VTKPolyhedron) = VTKCellTypes.VTK_POLYHEDRON.vtk_id % UInt8
cell_type_id(::MeshCell{<:PolyData.CellType}) = 0x00  # PolyData has no Types dataset

# Connectivity/offset/type arrays for one partition of an unstructured grid.
# Polyhedron face information is collected separately (empty when no
# polyhedra are present). All indices converted to 0-based.
struct PartitionTopology
    connectivity::Vector{Int64}
    offsets::Vector{Int64}       # ncells + 1, starting at 0
    types::Vector{UInt8}
    # polyhedron support
    n_polyhedra::Int
    face_connectivity::Vector{Int64}
    face_offsets::Vector{Int64}          # nfaces + 1, starting at 0
    polyhedron_to_faces::Vector{Int64}   # partition-local face ids
    polyhedron_offsets::Vector{Int64}    # ncells + 1, starting at 0
end

function PartitionTopology(cells::AbstractVector{<:AnyCell}, npoints::Integer = typemax(Int))
    nconn = sum(c -> length(c.connectivity), cells; init = 0)
    connectivity = Vector{Int64}(undef, nconn)
    offsets = Vector{Int64}(undef, length(cells) + 1)
    types = Vector{UInt8}(undef, length(cells))
    offsets[1] = 0
    n_polyhedra = 0
    face_connectivity = Int64[]
    face_offsets = Int64[0]
    polyhedron_to_faces = Int64[]
    polyhedron_offsets = Vector{Int64}(undef, length(cells) + 1)
    polyhedron_offsets[1] = 0
    pos = 0
    for (i, cell) in enumerate(cells)
        for id in cell.connectivity
            1 <= id <= npoints || throw(ArgumentError("cell $i references point $id, valid range is 1:$npoints"))
            pos += 1
            connectivity[pos] = id - 1
        end
        offsets[i + 1] = pos
        types[i] = cell_type_id(cell)
        polyhedron_offsets[i + 1] = polyhedron_offsets[i]
        if cell isa VTKPolyhedron
            n_polyhedra += 1
            for face in VTKBase.faces(cell)
                for id in face
                    1 <= id <= npoints || throw(ArgumentError("polyhedron $i face references point $id, valid range is 1:$npoints"))
                    push!(face_connectivity, id - 1)
                end
                push!(face_offsets, length(face_connectivity))
                push!(polyhedron_to_faces, length(face_offsets) - 2)  # 0-based face id
            end
            polyhedron_offsets[i + 1] += length(VTKBase.faces(cell))
        end
    end
    return PartitionTopology(
        connectivity, offsets, types,
        n_polyhedra, face_connectivity, face_offsets,
        polyhedron_to_faces, polyhedron_offsets
    )
end

n_cells(t::PartitionTopology) = length(t.types)
n_faces(t::PartitionTopology) = length(t.face_offsets) - 1

# Points input conversion: always a (3, N) Float matrix (2-D/1-D padded with
# zeros), written to a HDF (N, 3) dataset.
function prepare_points(points::AbstractMatrix{T}) where {T <: Real}
    dim, n = size(points)
    Tp = T <: AbstractFloat ? T : Float64
    dim == 3 && return (n, Tp === T ? points : Tp.(points))
    1 <= dim <= 3 || throw(ArgumentError("points matrix must be (dim ≤ 3) × N, got $(size(points))"))
    padded = zeros(Tp, 3, n)
    padded[1:dim, :] .= points
    return (n, padded)
end

function prepare_points(points::AbstractVector{T}) where {T}
    if T <: Real  # a single coordinate vector: 1-D points
        return prepare_points(reshape(points, 1, :))
    end
    ncomp, n, mat = prepare_tuples(points)
    return prepare_points(ncomp == 1 ? reshape(mat, 1, :) : mat)
end

function prepare_points(xyz::Tuple{Vararg{AbstractVector{<:Real}}})
    length(xyz) <= 3 || throw(ArgumentError("at most three coordinate vectors expected"))
    n = length(xyz[1])
    all(v -> length(v) == n, xyz) || throw(ArgumentError("coordinate vectors must have equal length"))
    mat = zeros(promote_type(map(eltype, xyz)...), 3, n)
    for (i, v) in enumerate(xyz)
        mat[i, :] .= v
    end
    return (n, mat)
end
