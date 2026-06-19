module Ellipsoids

using StaticArrays
using LinearAlgebra: norm

export FillableEllipsoid, FillableSphere, radii

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableEllipsoid{T, F, I} <: AbstractFillableShape

Axis-aligned ellipsoid centered at `center_xyz` with semi-axes `radii_xyz`.

The `fill_function` receives normalized local coordinates
`((x-cx)/rx, (y-cy)/ry, (z-cz)/rz)`, where the surface is at unit norm.

# Fields
- `center_xyz`: center in world space
- `radii_xyz`: semi-axis lengths along x, y, z
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableEllipsoid{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    radii_xyz::SVector{3, T}
    fill_function::F
    interpolation::I
end

"""
    FillableSphere(center, radius, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableEllipsoid`](@ref) with equal radii along all three axes.
"""
function FillableSphere(center::NTuple{3, T}, radius::T, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableEllipsoid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(_ -> radius, 3)), f, interpolation)
end

"""
    FillableEllipsoid(center, radii, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableEllipsoid`](@ref) with independent semi-axis lengths `radii`.
"""
function FillableEllipsoid(center::NTuple{3, T}, radii::NTuple{3, T}, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableEllipsoid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(radii), f, interpolation)
end

Types.center(e::FillableEllipsoid) = e.center_xyz

"""
    radii(shape) -> SVector{3}

Return the semi-axis lengths of an ellipsoid or ellipsoid-derived shape.
"""
radii(e::FillableEllipsoid) = e.radii_xyz
Types.interpolation(e::FillableEllipsoid) = e.interpolation

function Base.in(point, e::FillableEllipsoid{T}) where {T}
    unit = one(Base.promote_op(Base.:/, T, T))
    (((x - x₀) / r)^2 for (x, x₀, r) in zip(point, center(e), radii(e))) |> sum <= unit
end

function Base.fill(e::FillableEllipsoid{T}, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}) where {T}
    local_coords = ntuple(i -> (voxel_center_xyz[i] - e.center_xyz[i]) / e.radii_xyz[i], 3)
    return e.fill_function(local_coords)
end

# Approximate SDF (standard k0*(k0-1)/k1 formula). Not an exact distance bound.
function Types.sdf(e::FillableEllipsoid{T}, point::NTuple{3,T}) where {T}
    p = SVector{3,T}(point) - e.center_xyz
    r = e.radii_xyz
    k0 = norm(p ./ r)
    k1 = norm(p ./ (r .* r))
    return k0 * (k0 - one(T)) / k1
end

# k0*(k0-1)/k1 overestimates distance outside the shell, so has_exact_sdf is false.

end
