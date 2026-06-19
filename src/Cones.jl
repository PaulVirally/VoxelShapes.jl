module Cones

using StaticArrays

export FillableCone

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableCone{T, F, I} <: AbstractFillableShape

Axis-aligned cone or frustum centered at the midpoint between its two faces.

The radius varies linearly from `base_radius` at the negative end of the axis
(`-half_height`) to `top_radius` at the positive end (`+half_height`). Setting
`top_radius = 0` gives a true cone; unequal nonzero values give a frustum.

The longitudinal axis is selected by `axis` (1 = x, 2 = y, 3 = z).
The `fill_function` receives local coordinates `(radial_fraction, axial_fraction, 0)`,
where `radial_fraction = r / r_at_that_height` and
`axial_fraction = axial_offset / half_height ∈ [-1, 1]`.

No exact closed-form SDF is available; `has_exact_sdf` returns `false`.

# Fields
- `center_xyz`: center in world space
- `base_radius`: radius at the `−half_height` face
- `top_radius`: radius at the `+half_height` face
- `half_height`: half the length along the longitudinal axis
- `axis`: longitudinal axis index (1, 2, or 3)
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableCone{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    base_radius::T
    top_radius::T
    half_height::T
    axis::Int
    fill_function::F
    interpolation::I
end

"""
    FillableCone(center, base_radius, top_radius, half_height, fill_val; axis=3, interpolation=LinearInterpolation())

Construct a [`FillableCone`](@ref).
"""
function FillableCone(center::NTuple{3, T}, base_radius::T, top_radius::T, half_height::T, fill_val;
                      axis::Int=3, interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCone{T, typeof(f), I}(SVector{3,T}(center), base_radius, top_radius, half_height, axis, f, interpolation)
end

Types.center(c::FillableCone) = c.center_xyz
Types.interpolation(c::FillableCone) = c.interpolation

@inline function _cone_radial_axial(point::NTuple{3,T}, c::FillableCone{T}) where {T}
    ax = c.axis
    p = SVector{3,T}(point) - c.center_xyz
    axial = p[ax]
    i1 = ax == 1 ? 2 : 1
    i2 = ax == 3 ? 2 : 3
    radial = sqrt(p[i1]^2 + p[i2]^2)
    return radial, axial
end

@inline function _interp_radius(c::FillableCone{T}, axial::T) where {T}
    # axial in [-half_height, half_height]; t in [0,1]
    t = (axial + c.half_height) / (2 * c.half_height)
    return c.base_radius + t * (c.top_radius - c.base_radius)
end

function Base.in(point::NTuple{3,T}, c::FillableCone{T}) where {T}
    radial, axial = _cone_radial_axial(point, c)
    abs(axial) > c.half_height && return false
    return radial <= _interp_radius(c, axial)
end

function Base.fill(c::FillableCone{T}, voxel_center_xyz::NTuple{3,T}, voxel_size_xyz::NTuple{3,T}) where {T}
    radial, axial = _cone_radial_axial(voxel_center_xyz, c)
    r_at = _interp_radius(c, axial)
    radial_frac = r_at > zero(T) ? radial / r_at : zero(T)
    axial_frac = axial / c.half_height
    return c.fill_function((radial_frac, axial_frac, zero(T)))
end

# No exact closed-form SDF for a frustum.

end # module
