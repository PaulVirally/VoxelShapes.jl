module Capsules

using StaticArrays
using LinearAlgebra: dot, norm

export FillableCapsule

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableCapsule{T, F, I} <: AbstractFillableShape

Rounded rod (Minkowski sum of a line segment and a sphere) defined by two
endpoints `a`, `b` and a tube radius.

A point is inside when its distance to segment `ab` is at most `radius`.
The `fill_function` receives local coordinates `(dist_to_segment/radius, t, 0)`,
where `t ∈ [0, 1]` is the projection parameter along the segment.

Has an exact SDF.

# Fields
- `a`: first segment endpoint in world space
- `b`: second segment endpoint in world space
- `radius`: tube radius
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableCapsule{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    a::SVector{3, T}
    b::SVector{3, T}
    radius::T
    fill_function::F
    interpolation::I
end

"""
    FillableCapsule(a, b, radius, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableCapsule`](@ref) from segment endpoints `a` and `b`.
"""
function FillableCapsule(a::NTuple{3, T}, b::NTuple{3, T}, radius::T, fill_val;
                         interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCapsule{T, typeof(f), I}(SVector{3,T}(a), SVector{3,T}(b), radius, f, interpolation)
end

Types.interpolation(c::FillableCapsule) = c.interpolation

# Returns (distance_to_segment, t_parameter) where t ∈ [0,1] along ab.
@inline function _segment_dist(point::NTuple{3,T}, c::FillableCapsule{T}) where {T}
    p = SVector{3,T}(point)
    ab = c.b - c.a
    ap = p - c.a
    len2 = dot(ab, ab)
    t = len2 > zero(T) ? clamp(dot(ap, ab) / len2, zero(T), one(T)) : zero(T)
    closest = c.a + t * ab
    return norm(p - closest), t
end

function Base.in(point::NTuple{3,T}, c::FillableCapsule{T}) where {T}
    dist, _ = _segment_dist(point, c)
    return dist <= c.radius
end

function Base.fill(c::FillableCapsule{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    dist, t = _segment_dist(voxel_center_xyz, c)
    local_coords = (dist / c.radius, t, zero(T))
    return c.fill_function(local_coords)
end

# Exact SDF for a capsule
function Types.sdf(c::FillableCapsule{T}, point::NTuple{3,T}) where {T}
    dist, _ = _segment_dist(point, c)
    return dist - c.radius
end

Types.has_exact_sdf(::FillableCapsule) = true

end # module
