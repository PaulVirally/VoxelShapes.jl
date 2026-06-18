module Ellipsoids

using StaticArrays

export FillableEllipsoid, FillableSphere, center, radii

using ..Types
using ..Interpolations: LinearInterpolation

struct FillableEllipsoid{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    radii_xyz::SVector{3, T}
    fill_function::F
    interpolation::I
end

function FillableSphere(center::NTuple{3, T}, radius::T, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableEllipsoid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(_ -> radius, 3)), f, interpolation)
end

function FillableEllipsoid(center::NTuple{3, T}, radii::NTuple{3, T}, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableEllipsoid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(radii), f, interpolation)
end

center(e::FillableEllipsoid) = e.center_xyz
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

end
