module Slabs

using StaticArrays
using LinearAlgebra: dot, norm

export FillableSlab, FillableHalfSpace

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableHalfSpace{T, F, I} <: AbstractFillableShape

The half-space on the inward side of a plane.

A point `p` is inside when `dot(p - point, normal) ≤ 0`, where `normal` is the
outward unit normal. The `fill_function` receives local coordinates
`(signed_distance, 0, 0)` where `signed_distance = dot(p - point, normal)`.

Has an exact SDF.

# Fields
- `point`: any point on the bounding plane
- `normal`: outward unit normal (normalized on construction)
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableHalfSpace{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    point::SVector{3, T}
    normal::SVector{3, T}   # unit normal (outward)
    fill_function::F
    interpolation::I
end

"""
    FillableHalfSpace(point, normal, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableHalfSpace`](@ref). `normal` is normalized automatically.
"""
function FillableHalfSpace(point::NTuple{3, T}, normal::NTuple{3, T}, fill_val;
                           interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    n = SVector{3,T}(normal)
    n = n / norm(n)
    f = _ -> fill_val
    FillableHalfSpace{T, typeof(f), I}(SVector{3,T}(point), n, f, interpolation)
end

Types.interpolation(h::FillableHalfSpace) = h.interpolation

function Base.in(point::NTuple{3,T}, h::FillableHalfSpace{T}) where {T}
    dot(SVector{3,T}(point) - h.point, h.normal) <= zero(T)
end

function Base.fill(h::FillableHalfSpace{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    d = dot(SVector{3,T}(voxel_center_xyz) - h.point, h.normal)
    return h.fill_function((d, zero(T), zero(T)))
end

# Exact SDF: signed distance to the plane (positive outside = in direction of normal)
Types.sdf(h::FillableHalfSpace{T}, point::NTuple{3,T}) where {T} =
    dot(SVector{3,T}(point) - h.point, h.normal)

Types.has_exact_sdf(::FillableHalfSpace) = true

"""
    FillableSlab{T, F, I} <: AbstractFillableShape

Infinite planar slab: points within `half_thickness` of a reference plane.

A point `p` is inside when `|dot(p - point, normal)| ≤ half_thickness`. The
`fill_function` receives local coordinates `(d / half_thickness, 0, 0)`, where
`d = dot(p - point, normal)`, so the coordinate ranges from -1 to +1 across
the slab.

Has an exact SDF.

# Fields
- `point`: any point on the midplane
- `normal`: unit normal to the slab (normalized on construction)
- `half_thickness`: half the slab thickness
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableSlab{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    point::SVector{3, T}
    normal::SVector{3, T}   # unit normal
    half_thickness::T
    fill_function::F
    interpolation::I
end

"""
    FillableSlab(point, normal, half_thickness, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableSlab`](@ref). `normal` is normalized automatically.
"""
function FillableSlab(point::NTuple{3, T}, normal::NTuple{3, T}, half_thickness::T, fill_val;
                      interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    n = SVector{3,T}(normal)
    n = n / norm(n)
    f = _ -> fill_val
    FillableSlab{T, typeof(f), I}(SVector{3,T}(point), n, half_thickness, f, interpolation)
end

Types.interpolation(s::FillableSlab) = s.interpolation

function Base.in(point::NTuple{3,T}, s::FillableSlab{T}) where {T}
    abs(dot(SVector{3,T}(point) - s.point, s.normal)) <= s.half_thickness
end

function Base.fill(s::FillableSlab{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    d = dot(SVector{3,T}(voxel_center_xyz) - s.point, s.normal)
    return s.fill_function((d / s.half_thickness, zero(T), zero(T)))
end

# Exact SDF for a slab
Types.sdf(s::FillableSlab{T}, point::NTuple{3,T}) where {T} =
    abs(dot(SVector{3,T}(point) - s.point, s.normal)) - s.half_thickness

Types.has_exact_sdf(::FillableSlab) = true

end # module
