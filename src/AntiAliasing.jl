module AntiAliasing

export aa, NoAntiAliasing, SuperResolutionAntiAliasing, GaussianAntiAliasing
export SubpixelAntiAliasing, AdaptiveAntiAliasing

using StaticArrays
using ..Types
using ..Interpolations: interp_init, interp_accumulate, interp_finalize

"""
    aa(shape, voxel_center_xyz, voxel_size_xyz, background, anti_alias) -> value

Compute the rasterized value for a single voxel given an anti-aliasing strategy.

Returns `background` when the voxel is outside the shape, or a blended value
on the boundary according to the chosen `anti_alias` method.
"""
function aa end

"""
    NoAntiAliasing <: AbstractAntiAliasing

Point-sampled rasterization with no blending.

Each voxel is either fully inside the shape (returns `fill(shape, ...)`) or
fully outside (returns `background`), tested at the voxel center.
"""
struct NoAntiAliasing <: AbstractAntiAliasing end

function aa(shape::AbstractFillableShape, voxel_center_xyz::NTuple{3, T}, voxel_size_xyz::NTuple{3, T}, background::U, ::NoAntiAliasing) where {T, U}
    if voxel_center_xyz in shape
        return fill(shape, voxel_center_xyz, voxel_size_xyz)
    end
    return background
end

"""
    SuperResolutionAntiAliasing <: AbstractAntiAliasing

Anti-aliasing by averaging over a regular sub-grid within each voxel.

# Fields
- `super_grid`: number of sub-samples along each axis. A value of `(n, n, n)`
  produces n³ samples per voxel.

Construct with `SuperResolutionAntiAliasing(n)` for an isotropic n³ grid,
or `SuperResolutionAntiAliasing((nx, ny, nz))` for an anisotropic one.
"""
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

"""
    GaussianAntiAliasing{Nx, Ny, Nz, T} <: AbstractAntiAliasing

Anti-aliasing using a separable Gaussian kernel.

The kernel is pre-normalized so weights sum to one. Each voxel is sampled
at `Nx × Ny × Nz` points spaced one voxel apart, weighted by the product of
per-axis Gaussian weights.

Construct with `GaussianAntiAliasing(σ, kernel_size)` for isotropic smoothing,
or `GaussianAntiAliasing((σx, σy, σz), (Nx, Ny, Nz))` for anisotropic.
Kernel sizes must be odd.
"""
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

"""
    SubpixelAntiAliasing <: AbstractAntiAliasing

Analytic coverage estimate using the signed distance function. O(1) per voxel,
GPU-safe.

Approximates the fraction of the voxel inside the shape as a linear ramp through
the SDF value at the voxel center, normalized by half the voxel diagonal. The
boundary is assumed locally planar; accuracy degrades when voxel size approaches
the surface's radius of curvature. Coverage is slightly approximate for
anisotropic voxels.

Requires `sdf(shape, point)` to be implemented.
"""
struct SubpixelAntiAliasing <: AbstractAntiAliasing end

function aa(shape::AbstractFillableShape, vc::NTuple{3,T}, vs::NTuple{3,T}, bg::U, ::SubpixelAntiAliasing) where {T, U}
    d = sdf(shape, vc)
    h = oftype(d, 0.5) * sqrt(vs[1]^2 + vs[2]^2 + vs[3]^2)
    frac = clamp(oftype(d, 0.5) - d / (2 * h), zero(d), one(d))
    frac <= zero(d) && return bg
    interp = interpolation(shape)
    acc = interp_init(interp, U)
    acc = interp_accumulate(interp, acc, fill(shape, vc, vs), frac)
    acc = interp_accumulate(interp, acc, bg, one(frac) - frac)
    return interp_finalize(interp, acc)
end

"""
    AdaptiveAntiAliasing{A<:AbstractAntiAliasing} <: AbstractAntiAliasing

Wrapper that skips the inner anti-aliasing stencil for voxels clearly inside or
outside the shape.

When `has_exact_sdf(shape)` is `true`, a single SDF evaluation determines
whether the voxel is at least half a diagonal away from the surface. Only
boundary voxels fall through to the wrapped `inner` strategy.

# Fields
- `inner`: the anti-aliasing method applied to boundary voxels
"""
struct AdaptiveAntiAliasing{A<:AbstractAntiAliasing} <: AbstractAntiAliasing
    inner::A
end

function aa(shape::AbstractFillableShape, vc::NTuple{3,T}, vs::NTuple{3,T}, bg::U, a::AdaptiveAntiAliasing) where {T, U}
    if has_exact_sdf(shape)
        d = sdf(shape, vc)
        h = oftype(d, 0.5) * sqrt(vs[1]^2 + vs[2]^2 + vs[3]^2)
        d <= -h && return fill(shape, vc, vs)
        d >= h && return bg
    end
    return aa(shape, vc, vs, bg, a.inner)
end

end
