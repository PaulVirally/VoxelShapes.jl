module Types

export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation
export interpolation, sdf, has_exact_sdf, center

"""
    AbstractFillableShape

Base type for all shapes that can be rasterized into a [`World`](@ref VoxelShapes.World).

Subtypes must implement:
- `Base.in(point, shape)` (containment test)
- `Base.fill(shape, voxel_center, voxel_size)` (fill value at a voxel center)
- `interpolation(shape)` (the interpolation scheme used during anti-aliasing)
- `sdf(shape, point)` (signed distance to the surface (negative inside))
"""
abstract type AbstractFillableShape end

"""
    AbstractAntiAliasing

Base type for anti-aliasing strategies applied per voxel.

Concrete subtypes: `NoAntiAliasing`, `SuperResolutionAntiAliasing`, `GaussianAntiAliasing`,
`SubpixelAntiAliasing`, `AdaptiveAntiAliasing`.
"""
abstract type AbstractAntiAliasing end

"""
    AbstractInterpolation

Base type for accumulation strategies used when blending sub-voxel samples.

Subtypes implement the three-function fold protocol: `interp_init`, `interp_accumulate`, `interp_finalize`.
"""
abstract type AbstractInterpolation end

"""
    interpolation(shape) -> AbstractInterpolation

Return the interpolation scheme attached to `shape`.
"""
function interpolation end

"""
    center(shape) -> SVector{3}

Return the geometric center of `shape`.

Not all shapes define a meaningful center; those that do override this method.
"""
function center end

"""
    sdf(shape, point::NTuple{3,T}) -> T

Signed distance from `point` to the boundary of `shape`.

Negative inside the shape, positive outside, zero on the surface.
Shapes where this is an exact Euclidean distance report `has_exact_sdf(shape) == true`.
"""
function sdf end

"""
    has_exact_sdf(shape) -> Bool

Return `true` if `sdf(shape, ...)` is a true Euclidean signed distance function.

The default is `false`. Shapes with exact SDFs allow `AdaptiveAntiAliasing`
to skip the inner stencil for voxels clearly inside or outside the surface.
"""
has_exact_sdf(::AbstractFillableShape) = false

end
