# Extending

A custom shape needs four methods:

```julia
# 1. Containment test
Base.in(point::NTuple{3,T}, shape::MyShape) -> Bool

# 2. Fill value at a voxel center
Base.fill(shape::MyShape, voxel_center::NTuple{3,T}, voxel_size::NTuple{3,T}) -> value

# 3. Interpolation strategy
VoxelShapes.interpolation(shape::MyShape) -> AbstractInterpolation

# 4. Signed distance function
VoxelShapes.sdf(shape::MyShape, point::NTuple{3,T}) -> T
```

If your SDF is a true Euclidean distance (negative inside, zero on the
surface, positive outside), also add:

```julia
VoxelShapes.has_exact_sdf(::MyShape) = true
```

This lets [`AdaptiveAntiAliasing`](@ref) skip the inner stencil for voxels
far from the surface.

Optionally, add a bounding box:

```julia
VoxelShapes.bounding_box(shape::MyShape) -> (lower::NTuple{3,T}, upper::NTuple{3,T})
```

The default returns the all-infinite box `((-Inf,-Inf,-Inf),
(Inf,Inf,Inf))`. [`refine`](@ref) around a shape without this method
refines the whole grid. A tighter, conservative box (never smaller than
the shape, larger is fine) lets `refine` carve out just the region around
the shape. See [Grids and refinement](@ref) for how the box is used.

## Minimal example

```julia
using VoxelShapes
using VoxelShapes: AbstractFillableShape, AbstractInterpolation

struct MyBox{T, I<:AbstractInterpolation} <: AbstractFillableShape
    center::NTuple{3,T}
    half_size::T
    fill_val::T
    interpolation::I
end

Base.in(p::NTuple{3,T}, s::MyBox{T}) where T =
    all(abs(p[i] - s.center[i]) <= s.half_size for i in 1:3)

Base.fill(s::MyBox{T}, vc::NTuple{3,T}, vs::NTuple{3,T}) where T = s.fill_val

VoxelShapes.interpolation(s::MyBox) = s.interpolation

function VoxelShapes.sdf(s::MyBox{T}, p::NTuple{3,T}) where T
    # Exact box SDF
    q = ntuple(i -> abs(p[i] - s.center[i]) - s.half_size, 3)
    pos = ntuple(i -> max(q[i], zero(T)), 3)
    return sqrt(sum(x^2 for x in pos)) + min(maximum(q), zero(T))
end

VoxelShapes.has_exact_sdf(::MyBox) = true
```

## GPU compatibility

For a shape to work with `rasterize(geometry, region, CuArray)`, the
struct must be `isbits`. Avoid heap-allocated fields: arrays, strings,
arbitrary closures. Use `StaticArrays.SVector` for fixed-size vectors and
`NTuple` for fixed-size tuples.

See the [API reference](@ref "API reference") for the abstract type documentation.
