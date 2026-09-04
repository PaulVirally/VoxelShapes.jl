module Cylinders

using StaticArrays
using LinearAlgebra: norm

export FillableCylinder, radius, half_height

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableCylinder{T, F, I} <: AbstractFillableShape

Finite axis-aligned cylinder centered at `center_xyz`.

# Fields
- `center_xyz`: center in world space
- `radius`: circular cross-section radius
- `half_height`: half the length along the longitudinal axis
- `axis`: longitudinal axis index (1, 2, or 3)
- `fill_function`: maps local coords `(radial_fraction, axial_fraction, 0)`
  to a fill value, both ranges `∈ [0, 1]` and `[-1, 1]`
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableCylinder{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    radius::T
    half_height::T
    axis::Int
    fill_function::F
    interpolation::I
end

"""
    FillableCylinder(center, radius, half_height, fill_val; axis=3, interpolation=LinearInterpolation())

Construct a [`FillableCylinder`](@ref).
"""
function FillableCylinder(center::NTuple{3, T}, radius::T, half_height::T, fill_val;
                          axis::Int=3, interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCylinder{T, typeof(f), I}(SVector{3,T}(center), radius, half_height, axis, f, interpolation)
end

Types.center(c::FillableCylinder) = c.center_xyz

"""
    radius(shape) -> T

Return the circular cross-section radius of a cylinder or capsule.
"""
radius(c::FillableCylinder) = c.radius

"""
    half_height(shape) -> T

Return half the length of a cylinder or cone along its longitudinal axis.
"""
half_height(c::FillableCylinder) = c.half_height
Types.interpolation(c::FillableCylinder) = c.interpolation

@inline function _radial_axial(point::NTuple{3,T}, c::FillableCylinder{T}) where {T}
    ax = c.axis
    p = SVector{3,T}(point) - c.center_xyz
    axial = p[ax]
    # radial distance in the two non-axis dims
    i1 = ax == 1 ? 2 : 1
    i2 = ax == 3 ? 2 : 3
    radial = sqrt(p[i1]^2 + p[i2]^2)
    return radial, axial
end

function Base.in(point::NTuple{3,T}, c::FillableCylinder{T}) where {T}
    radial, axial = _radial_axial(point, c)
    return radial <= c.radius && abs(axial) <= c.half_height
end

function Base.fill(c::FillableCylinder{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    radial, axial = _radial_axial(voxel_center_xyz, c)
    local_coords = (radial / c.radius, axial / c.half_height, zero(T))
    return c.fill_function(local_coords)
end

# Exact SDF for an axis-aligned finite cylinder
function Types.sdf(c::FillableCylinder{T}, point::NTuple{3,T}) where {T}
    radial, axial = _radial_axial(point, c)
    d = SVector(radial - c.radius, abs(axial) - c.half_height)
    return min(maximum(d), zero(T)) + norm(max.(d, zero(T)))
end

Types.has_exact_sdf(::FillableCylinder) = true

"""
    bounding_box(shape::FillableCylinder) -> (lower, upper)

Exact axis-aligned bounding box: `± half_height` along `axis`, `± radius` on
the other two axes.
"""
function Types.bounding_box(c::FillableCylinder{T}) where {T}
    ax = c.axis
    ctr = c.center_xyz
    lower = ntuple(i -> ctr[i] - (i == ax ? c.half_height : c.radius), 3)
    upper = ntuple(i -> ctr[i] + (i == ax ? c.half_height : c.radius), 3)
    return (lower, upper)
end

end # module
