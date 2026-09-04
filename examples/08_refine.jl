# 08_refine.jl
#
# A CompositeGrid glues Regions of different resolutions into one grid. This
# example refines a coarse grid around a small sphere and rasterizes into
# the composite grid. It then looks at the result two ways: as a uniform
# array (regrid, for plotting), and colored by region index, to show the
# tiling itself.

using VoxelShapes
using GLMakie

N = 16
coarse = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))

sphere = FillableSphere((0.5, 0.5, 0.5), 0.12, 1.0)

# Refine by a factor of 2 around the sphere's bounding box. Add one coarse
# voxel of padding so the anti-aliased boundary voxels land in the fine core.
grid = refine(coarse, sphere; factor=2, padding=1 // N)

geometry = Geometry([sphere], 0.0, AdaptiveAntiAliasing(SuperResolutionAntiAliasing(4)))
field = rasterize(geometry, grid)

# A dense array at the finest scale present in the grid, for plotting.
uniform = regrid(field)

# A field of region indices, not rasterized values, regridded the same way.
# Each block of uniform color is one region, showing the tiling itself.
index_field = CompositeField(
    reduce(vcat, fill(Float64(idx), length(reg)) for (idx, reg) in enumerate(regions(grid))),
    grid,
)
region_map = regrid(index_field)

fig = Figure(size = (900, 420))

nz = size(uniform, 3)
ax1 = Axis(fig[1, 1], title = "Sphere, refined core, center slice", aspect = DataAspect())
heatmap!(ax1, uniform[:, :, nz ÷ 2], colormap = :grays, colorrange = (0, 1))
hidedecorations!(ax1)

ax2 = Axis(fig[1, 2], title = "Region index, same slice", aspect = DataAspect())
heatmap!(ax2, region_map[:, :, nz ÷ 2], colormap = :tab10)
hidedecorations!(ax2)

save(joinpath(@__DIR__, "08_refine.png"), fig)
display(fig)
