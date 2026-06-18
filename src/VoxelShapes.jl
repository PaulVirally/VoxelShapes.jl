module VoxelShapes

export World

using KernelAbstractions
using CUDA

include("Types.jl")
using .Types
export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation, interpolation

include("Interpolations.jl")
using .Interpolations
export LinearInterpolation, interp_init, interp_accumulate, interp_finalize

include("Ellipsoids.jl")
using .Ellipsoids
export FillableEllipsoid, FillableSphere, center, radii

include("Cuboids.jl")
using .Cuboids
export FillableCuboid, FillableCube, center, half_lengths, lengths

include("AntiAliasing.jl")
using .AntiAliasing
export aa, NoAntiAliasing, SuperResolutionAntiAliasing, GaussianAntiAliasing

struct World{T, U}
    num_voxels_xyz::NTuple{3, Int}
    voxel_size_xyz::NTuple{3, T} # Convention is (0, 0, 0) is the front-down-left corner of the world
    shapes::AbstractVector{<:AbstractFillableShape}
    background::U # Default value where no shape fills the voxel
    aa::AbstractAntiAliasing
end

Base.size(world::World) = world.num_voxels_xyz

idx2global_xyz(idx::CartesianIndex, world::World) = Tuple((idx[i] - 1) * world.voxel_size_xyz[i] for i in 1:3)

@kernel function fill_voxel!(arr::AbstractArray{T, 3}, world::World) where {T}
    idx = @index(Global, Cartesian)
    arr[idx] = world.background # Default value if no shape fills the voxel
    fdl_corner_xyz = idx2global_xyz(idx, world) # Global coordinates of the voxel's front-down-left corner
    half = one(eltype(world.voxel_size_xyz)) / 2
    voxel_center_xyz = Tuple(fdl_corner_xyz[i] + half * world.voxel_size_xyz[i] for i in 1:3)
    for shape in world.shapes
        fill_value = aa(shape, voxel_center_xyz, world.voxel_size_xyz, world.background, world.aa)
        if fill_value !== world.background
            # If the shape fills the voxel, set the value in the array
            arr[idx] = fill_value
            break
        end
    end
end

function Base.fill!(arr::Array{T, 3}, world::World) where {T}
    backend = CPU()
    kernel = fill_voxel!(backend)
    kernel(arr, world, ndrange=size(arr))
end

function Base.fill!(arr::CuArray{T, 3}, world::World) where {T}
    backend = CUDABackend()
    if !(world.shapes isa Vector)
        shapes = CuArray(world.shapes)
    end
    kernel = fill_voxel!(backend)
    kernel(arr, world, ndrange=size(arr))
end

Base.push!(world::World, shape::AbstractFillableShape) = push!(world.shapes, shape)

function Base.Array(world::World)
    arr = zeros(eltype(world.background), world.num_voxels_xyz)
    fill!(arr, world)
    return arr
end

function CUDA.CuArray(world::World)
    arr = CUDA.zeros(eltype(world.background), world.num_voxels_xyz)
    fill!(arr, world)
    return arr
end

end # module VoxelShapes
