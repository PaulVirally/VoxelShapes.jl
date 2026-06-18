module AntiAliasing

export aa, NoAntiAliasing, SuperResolutionAntiAliasing, GaussianAntiAliasing

using StaticArrays
using ..Types
using ..Interpolations: interp_init, interp_accumulate, interp_finalize

struct NoAntiAliasing <: AbstractAntiAliasing end

function aa(shape::AbstractFillableShape, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}, background::U, ::NoAntiAliasing) where {T, U}
    if voxel_center_xyz in shape
        return fill(shape, voxel_center_xyz, voxel_size_xyz)
    end
    return background
end

struct SuperResolutionAntiAliasing <: AbstractAntiAliasing
    super_grid::NTuple{3, Int}
end

SuperResolutionAntiAliasing(super_resolution::Int=10) = SuperResolutionAntiAliasing((ntuple(_ -> super_resolution, 3)))

function aa(shape::AbstractFillableShape, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}, background::U, anti_alias::SuperResolutionAntiAliasing) where {T, U}
    interp = interpolation(shape)
    acc = interp_init(interp, U)
    grid = anti_alias.super_grid
    weight = one(U) / prod(grid)
    super_grid_stride = Tuple(voxel_size_xyz[i] / T(grid[i]) for i in 1:3)

    half = one(T) / 2
    offset = Tuple(-half * (one(T) - one(T) / grid[i]) * voxel_size_xyz[i] for i in 1:3)

    for i in 1:grid[1], j in 1:grid[2], k in 1:grid[3]
        super_voxel_center_xyz = (
            voxel_center_xyz[1] + offset[1] + (i - half) * super_grid_stride[1],
            voxel_center_xyz[2] + offset[2] + (j - half) * super_grid_stride[2],
            voxel_center_xyz[3] + offset[3] + (k - half) * super_grid_stride[3]
        )
        val = aa(shape, Tuple(super_voxel_center_xyz), super_grid_stride, background, NoAntiAliasing())
        acc = interp_accumulate(interp, acc, val, weight)
    end
    return interp_finalize(interp, acc)
end

struct GaussianAntiAliasing{Nx, Ny, Nz, T} <: AbstractAntiAliasing
    kernel_data::Tuple{SVector{Nx, T}, SVector{Ny, T}, SVector{Nz, T}}
end

function GaussianAntiAliasing(σs::NTuple{3, T}, kernel_size::NTuple{3, Int}) where {T}
    @assert all(isodd, kernel_size) "Kernel size must be odd in all dimensions"
    half_kernel_sizes = Tuple(kernel_size[i] ÷ 2 for i in 1:3)

    kernels = (SVector{kernel_size[i], T}(exp(-T(x^2) / (2σs[i]^2)) for x in -(half_kernel_sizes[i]):(half_kernel_sizes[i])) for i in 1:3)
    kernels = (k / sum(k) for k in kernels)

    return GaussianAntiAliasing{kernel_size[1], kernel_size[2], kernel_size[3], T}(Tuple(kernels))
end
GaussianAntiAliasing(σ::T, kernel_size::Int) where {T} = GaussianAntiAliasing((σ, σ, σ), (kernel_size, kernel_size, kernel_size))

function aa(shape::AbstractFillableShape, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}, background::U, anti_alias::GaussianAntiAliasing{Nx, Ny, Nz, U}) where {T, U, Nx, Ny, Nz}
    interp = interpolation(shape)
    acc = interp_init(interp, U)
    kernel_x, kernel_y, kernel_z = anti_alias.kernel_data
    half_size_x, half_size_y, half_size_z = (Nx ÷ 2, Ny ÷ 2, Nz ÷ 2)
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        dx, dy, dz = i - half_size_x - 1, j - half_size_y - 1, k - half_size_z - 1
        super_voxel_center_xyz = (
            voxel_center_xyz[1] + T(dx) * voxel_size_xyz[1] / T(Nx),
            voxel_center_xyz[2] + T(dy) * voxel_size_xyz[2] / T(Ny),
            voxel_center_xyz[3] + T(dz) * voxel_size_xyz[3] / T(Nz)
        )
        val = aa(shape, super_voxel_center_xyz, voxel_size_xyz, background, NoAntiAliasing())
        weight = kernel_x[i] * kernel_y[j] * kernel_z[k]
        acc = interp_accumulate(interp, acc, val, weight)
    end
    return interp_finalize(interp, acc)
end

end
