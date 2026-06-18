module Cuboids

using StaticArrays

export FillableCuboid, FillableCube, center, half_lengths, lengths

using ..Types
using ..Interpolations: LinearInterpolation

struct FillableCuboid{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    half_lengths_xyz::SVector{3, T}
    fill_function::F
    interpolation::I
end

function FillableCube(center::NTuple{3, T}, length::T, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCuboid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(_ -> length/2, 3)), f, interpolation)
end

function FillableCuboid(center::NTuple{3, T}, lengths::NTuple{3, T}, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCuboid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(i -> lengths[i]/2, 3)), f, interpolation)
end

center(c::FillableCuboid) = c.center_xyz
half_lengths(c::FillableCuboid) = c.half_lengths_xyz
lengths(c::FillableCuboid) = 2 .* c.half_lengths_xyz
Types.interpolation(c::FillableCuboid) = c.interpolation

function Base.in(point::NTuple{3,T}, c::FillableCuboid{T}) where {T}
    all(abs(x - x₀) <= hl for (x, x₀, hl) in zip(point, center(c), half_lengths(c)))
end

function Base.fill(c::FillableCuboid{T}, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}) where {T}
    local_coords = ntuple(i -> (voxel_center_xyz[i] - c.center_xyz[i]) / c.half_lengths_xyz[i], 3)
    return c.fill_function(local_coords)
end

end # module
