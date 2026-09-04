module Types

export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation
export interpolation, sdf, has_exact_sdf, center, bounding_box

"""
    AbstractFillableShape

Base type for all shapes that can be placed in a [`Geometry`](@ref VoxelShapes.Geometry) and rasterized.

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

"""
    bounding_box(shape) -> (lower::NTuple{3,T}, upper::NTuple{3,T})

Return an axis-aligned bounding box for `shape` as a `(lower, upper)` pair of
`NTuple{3,T}` corners.

The box is conservative: it may be larger than `shape`, but never smaller
(every point in `shape` lies within it). Infinite extents are reported as
`±Inf`.

The default method, used for shapes that do not override it, returns the
all-infinite box `((-Inf,-Inf,-Inf), (Inf,Inf,Inf))`. Refining a grid around
such a shape then refines everything rather than failing.
"""
bounding_box(::AbstractFillableShape) = ((-Inf, -Inf, -Inf), (Inf, Inf, Inf))

end
