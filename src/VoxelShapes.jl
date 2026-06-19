module VoxelShapes

export World, add_shape

using KernelAbstractions
using CUDA

include("Types.jl")
using .Types
export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation
export interpolation, sdf, has_exact_sdf, center

include("Interpolations.jl")
using .Interpolations
export LinearInterpolation, HarmonicInterpolation, GeometricMeanInterpolation, MaxInterpolation, MinInterpolation, DielectricInterpolation, MetalInterpolation
export interp_init, interp_accumulate, interp_finalize

include("Ellipsoids.jl")
using .Ellipsoids
export FillableEllipsoid, FillableSphere, radii

include("Cuboids.jl")
using .Cuboids
export FillableCuboid, FillableCube, half_lengths, lengths

include("Cylinders.jl")
using .Cylinders
export FillableCylinder, radius, half_height

include("Slabs.jl")
using .Slabs
export FillableSlab, FillableHalfSpace

include("Cones.jl")
using .Cones
export FillableCone

include("Tori.jl")
using .Tori
export FillableTorus, major_radius, minor_radius

include("Capsules.jl")
using .Capsules
export FillableCapsule

include("Transforms.jl")
using .Transforms
export Rotated

include("CSG.jl")
using .CSG
export UnionShape, IntersectionShape, DifferenceShape, ComplementShape
export csg_union, csg_intersect, csg_diff, csg_complement

include("Fills.jl")
using .Fills
export ConstantFill, RadialGradient, AxialGradient

include("AntiAliasing.jl")
using .AntiAliasing
export aa, NoAntiAliasing, SuperResolutionAntiAliasing, GaussianAntiAliasing
export SubpixelAntiAliasing, AdaptiveAntiAliasing

"""
    World{T, U, S, A}

Immutable scene descriptor that holds the voxel grid geometry, a tuple of shapes,
a background value, and an anti-aliasing strategy.

Shapes are stored as a heterogeneous `Tuple` (not a `Vector`) so that
`Base.fill!` compiles a type-stable, GPU-safe loop over them. Shapes are
evaluated in order; the first one whose containment test succeeds claims the
voxel.

The `background` value also acts as a transparency sentinel: a shape that
covers a voxel but evaluates to a value equal to `background` is treated as
transparent there, and evaluation falls through to the next shape (or the
`background`). Filling with the background value therefore punches a hole
rather than painting that value over a lower layer.

Convert to an array with `Array(world)` (CPU) or `CuArray(world)` (CUDA).
Use [`add_shape`](@ref) to append a shape (returns a new `World`).

# Fields
- `num_voxels_xyz`: grid dimensions as `(nx, ny, nz)`
- `voxel_size_xyz`: physical size of one voxel along each axis
- `shapes`: tuple of `AbstractFillableShape` values, evaluated first-to-last
- `background`: fill value for voxels not covered by any shape
- `aa`: anti-aliasing strategy applied to every voxel
- `origin_xyz`: world-space position of the front-bottom-left corner of voxel (1,1,1)
"""
struct World{T, U, S<:Tuple, A<:AbstractAntiAliasing}
    num_voxels_xyz::NTuple{3, Int}
    voxel_size_xyz::NTuple{3, T}
    shapes::S
    background::U
    aa::A
    origin_xyz::NTuple{3, T}
end

# 5-arg constructor with default origin=(0,0,0). The auto-generated outer constructor
# handles the 6-arg Tuple case, so we only add the 5-arg shorthand here.
World(nv::NTuple{3,Int}, vs::NTuple{3,T}, shapes::S, bg::U, a::A) where {T,U,S<:Tuple,A<:AbstractAntiAliasing} =
    World{T,U,S,A}(nv, vs, shapes, bg, a, (zero(T),zero(T),zero(T)))

# Convenience constructor: accept a Vector and convert to Tuple
World(nv::NTuple{3,Int}, vs::NTuple{3,T}, shapes::AbstractVector, bg::U, a::A,
      origin::NTuple{3,T}=(zero(T),zero(T),zero(T))) where {T,U,A<:AbstractAntiAliasing} =
    World(nv, vs, Tuple(shapes), bg, a, origin)

Base.size(world::World) = world.num_voxels_xyz
Base.eltype(world::World) = typeof(world.background)

Base.similar(world::World) = zeros(eltype(world), world.num_voxels_xyz)

idx2global_xyz(idx::CartesianIndex, world::World) =
    ntuple(i -> world.origin_xyz[i] + (idx[i] - 1) * world.voxel_size_xyz[i], 3)

# Type-stable recursive fold over a Tuple of shapes (GPU-safe)
@inline _fill_shapes(::Tuple{}, vc, vs, bg, a) = bg
@inline function _fill_shapes(shapes::Tuple, vc, vs, bg, a)
    v = aa(first(shapes), vc, vs, bg, a)
    return v !== bg ? v : _fill_shapes(Base.tail(shapes), vc, vs, bg, a)
end

@kernel function fill_voxel!(arr::AbstractArray{T, 3}, world::World) where {T}
    idx = @index(Global, Cartesian)
    fdl_corner_xyz = idx2global_xyz(idx, world)
    half = one(eltype(world.voxel_size_xyz)) / 2
    voxel_center_xyz = ntuple(i -> fdl_corner_xyz[i] + half * world.voxel_size_xyz[i], 3)
    arr[idx] = _fill_shapes(world.shapes, voxel_center_xyz, world.voxel_size_xyz, world.background, world.aa)
end

function Base.fill!(arr::Array{T, 3}, world::World) where {T}
    kernel = fill_voxel!(CPU())
    kernel(arr, world, ndrange=size(arr))
end

function Base.fill!(arr::CuArray{T, 3}, world::World) where {T}
    kernel = fill_voxel!(CUDABackend())
    kernel(arr, world, ndrange=size(arr))
end

"""
    add_shape(world, shape) -> World

Return a new `World` with `shape` appended to the end of the shape tuple.

Because `World` is immutable and shapes are stored as a `Tuple`, this is the
functional replacement for `push!`.
"""
add_shape(world::World{T,U,S,A}, shape) where {T,U,S,A} =
    World(world.num_voxels_xyz, world.voxel_size_xyz,
          (world.shapes..., shape), world.background, world.aa, world.origin_xyz)

function Base.Array(world::World)
    arr = similar(world)
    fill!(arr, world)
    return arr
end

function CUDA.CuArray(world::World)
    arr = CUDA.zeros(eltype(world.background), world.num_voxels_xyz)
    fill!(arr, world)
    return arr
end

end # module VoxelShapes
