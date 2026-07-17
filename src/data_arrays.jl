# Conversion of user-provided data arrays to the (ncomponents, ntuples) form
# written to file (HDF5 row-major flips this to the (N, ncomp) layout VTK
# expects), plus data-location resolution.

# Returns (ncomp, ntuples, array). `array` is either a Vector (ncomp == 1,
# written as a 1-D dataset) or a (ncomp, ntuples) matrix.
prepare_tuples(A::AbstractVector{<:Real}) = (1, length(A), A)
prepare_tuples(A::AbstractMatrix{<:Real}) = (size(A, 1), size(A, 2), A)

# Vectors of isbits "tuple-like" element types (SVector, Tensors.Vec, NTuple, ...)
# with homogeneous field types are reinterpreted to a component matrix.
function prepare_tuples(A::AbstractVector{T}) where {T}
    isbitstype(T) || throw(ArgumentError("unsupported data eltype $T"))
    E = tuple_scalar_type(T)
    n = sizeof(T) ÷ sizeof(E)
    n * sizeof(E) == sizeof(T) || throw(ArgumentError("unsupported data eltype $T"))
    M = reinterpret(reshape, E, A)
    return n == 1 ? (1, length(A), vec(M)) : (n, length(A), M)
end

prepare_tuples(A::AbstractArray) =
    throw(ArgumentError("unsupported data array with ndims=$(ndims(A)); pass a vector or a (ncomponents, ntuples) matrix"))

function tuple_scalar_type(::Type{T}) where {T}
    T <: Real && return T
    fieldcount(T) == 0 && throw(ArgumentError("unsupported data eltype $T"))
    E = tuple_scalar_type(fieldtype(T, 1))
    for i in 2:fieldcount(T)
        tuple_scalar_type(fieldtype(T, i)) === E || throw(ArgumentError("unsupported data eltype $T (non-homogeneous fields)"))
    end
    return E
end

# Active-attribute marking (spec ≥ 1.0 group attributes + 2.6 array attribute).
const ATTRIBUTE_KINDS = (
    :Scalars, :Vectors, :Normals, :TCoords, :Tensors, :GlobalIds,
    :PedigreeIds, :EdgeFlag, :Tangents, :RationalWeights,
    :HigherOrderDegrees, :ProcessIds,
)

function mark_attribute(vtk, group, name::AbstractString, attribute::Symbol)
    attribute in ATTRIBUTE_KINDS ||
        throw(ArgumentError("unknown attribute kind $attribute; expected one of $(ATTRIBUTE_KINDS)"))
    write_ascii_attribute(group, String(attribute), String(name))
    write_ascii_attribute(group[name], "Attribute", String(attribute))
    bump_version!(vtk, (2, 6))
    return nothing
end
