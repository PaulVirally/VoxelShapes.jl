module Tori

using StaticArrays
using LinearAlgebra: norm

export FillableTorus, major_radius, minor_radius

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableTorus{T, F, I} <: AbstractFillableShape

Torus centered at `center_xyz`, with the ring lying in the plane perpendicular to `axis`.

`major_radius` (R) is the distance from the torus center to the tube center.
`minor_radius` (r) is the tube radius. A point is inside when
`(ρ - R)² + z² ≤ r²`, where ρ is the in-plane radial distance and z is the
axial distance from the center.

The longitudinal axis is selected by `axis` (1 = x, 2 = y, 3 = z; default 3
puts the ring in the xy-plane). The `fill_function` receives local coordinates
`(ρ/R, dist_to_tube_center/r, 0)`.

Has an exact SDF.

# Fields
- `center_xyz`: center in world space
- `major_radius`: ring radius (center to tube center)
- `minor_radius`: tube radius
- `axis`: symmetry axis index (1, 2, or 3)
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableTorus{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    major_radius::T
    minor_radius::T
    axis::Int
    fill_function::F
    interpolation::I
end

"""
    FillableTorus(center, major_radius, minor_radius, fill_val; axis=3, interpolation=LinearInterpolation())

Construct a [`FillableTorus`](@ref).
"""
function FillableTorus(center::NTuple{3, T}, major_radius::T, minor_radius::T, fill_val;
                       axis::Int=3, interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableTorus{T, typeof(f), I}(SVector{3,T}(center), major_radius, minor_radius, axis, f, interpolation)
end

Types.center(t::FillableTorus) = t.center_xyz

"""
    major_radius(shape) -> T

Return the major radius of a torus (distance from the torus center to the tube center).
"""
major_radius(t::FillableTorus) = t.major_radius

"""
    minor_radius(shape) -> T

Return the minor radius of a torus (tube radius).
"""
minor_radius(t::FillableTorus) = t.minor_radius
Types.interpolation(t::FillableTorus) = t.interpolation

@inline function _torus_components(point::NTuple{3,T}, t::FillableTorus{T}) where {T}
    ax = t.axis
    p = SVector{3,T}(point) - t.center_xyz
    axial = p[ax]
    i1 = ax == 1 ? 2 : 1
    i2 = ax == 3 ? 2 : 3
    rho = sqrt(p[i1]^2 + p[i2]^2)
    return rho, axial
end

function Base.in(point::NTuple{3,T}, t::FillableTorus{T}) where {T}
    rho, axial = _torus_components(point, t)
    return (rho - t.major_radius)^2 + axial^2 <= t.minor_radius^2
end

function Base.fill(t::FillableTorus{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    rho, axial = _torus_components(voxel_center_xyz, t)
    tube_dist = sqrt((rho - t.major_radius)^2 + axial^2)
    local_coords = (rho / t.major_radius, tube_dist / t.minor_radius, zero(T))
    return t.fill_function(local_coords)
end

# Exact SDF for a torus
function Types.sdf(t::FillableTorus{T}, point::NTuple{3,T}) where {T}
    rho, axial = _torus_components(point, t)
    return sqrt((rho - t.major_radius)^2 + axial^2) - t.minor_radius
end

Types.has_exact_sdf(::FillableTorus) = true

end # module
