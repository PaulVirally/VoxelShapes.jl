module Cuboids

using StaticArrays
using LinearAlgebra: norm

export FillableCuboid, FillableCube, half_lengths, lengths

using ..Types
using ..Interpolations: LinearInterpolation

"""
    FillableCuboid{T, F, I} <: AbstractFillableShape

Axis-aligned box centered at `center_xyz` with half-lengths `half_lengths_xyz`.

The `fill_function` receives local coordinates
`((x-cx)/hlx, (y-cy)/hly, (z-cz)/hlz)`, which range from -1 to +1 across the
box extent along each axis.

Has an exact SDF (Inigo Quilez formula).

# Fields
- `center_xyz`: center in world space
- `half_lengths_xyz`: half the box extent along each axis
- `fill_function`: callable mapping local coordinates to a fill value
- `interpolation`: blending strategy for anti-aliasing
"""
struct FillableCuboid{T, F, I<:AbstractInterpolation} <: AbstractFillableShape
    center_xyz::SVector{3, T}
    half_lengths_xyz::SVector{3, T}
    fill_function::F
    interpolation::I
end

"""
    FillableCube(center, length, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableCuboid`](@ref) with equal side lengths.
"""
function FillableCube(center::NTuple{3, T}, length::T, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCuboid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(_ -> length/2, 3)), f, interpolation)
end

"""
    FillableCuboid(center, lengths, fill_val; interpolation=LinearInterpolation())

Construct a [`FillableCuboid`](@ref) with independent side lengths `lengths`.
"""
function FillableCuboid(center::NTuple{3, T}, lengths::NTuple{3, T}, fill_val; interpolation::I=LinearInterpolation()) where {T, I<:AbstractInterpolation}
    f = _ -> fill_val
    FillableCuboid{T, typeof(f), I}(SVector{3,T}(center), SVector{3,T}(ntuple(i -> lengths[i]/2, 3)), f, interpolation)
end

Types.center(c::FillableCuboid) = c.center_xyz

"""
    half_lengths(shape) -> SVector{3}

Return the half-lengths (half the side lengths) of a cuboid.
"""
half_lengths(c::FillableCuboid) = c.half_lengths_xyz

"""
    lengths(shape) -> SVector{3}

Return the full side lengths of a cuboid.
"""
lengths(c::FillableCuboid) = 2 .* c.half_lengths_xyz
Types.interpolation(c::FillableCuboid) = c.interpolation

function Base.in(point::NTuple{3,T}, c::FillableCuboid{T}) where {T}
    all(abs(x - x₀) <= hl for (x, x₀, hl) in zip(point, center(c), half_lengths(c)))
end

function Base.fill(c::FillableCuboid{T}, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}) where {T}
    local_coords = ntuple(i -> (voxel_center_xyz[i] - c.center_xyz[i]) / c.half_lengths_xyz[i], 3)
    return c.fill_function(local_coords)
end

# Exact box SDF (Inigo Quilez formula)
function Types.sdf(c::FillableCuboid{T}, point::NTuple{3,T}) where {T}
    p = SVector{3,T}(point) - c.center_xyz
    q = abs.(p) - c.half_lengths_xyz
    return norm(max.(q, zero(T))) + min(maximum(q), zero(T))
end

Types.has_exact_sdf(::FillableCuboid) = true

end # module
