# 08_refine.jl
#
# A CompositeGrid glues Regions of different resolution into one grid. This
# refines a coarse grid around a small sphere, rasterizes into the composite
# grid, and looks at the result two ways: regridded to a uniform array (what
# you'd hand to a plotting function that expects a dense grid), and colored by
# region index instead, which shows the tiling itself rather than the field.

using VoxelShapes
using GLMakie

N = 16
coarse = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))

sphere = FillableSphere((0.5, 0.5, 0.5), 0.12, 1.0)

# Refine by a factor of 2 around the sphere's bounding box, with one coarse
# voxel of padding so the anti-aliased boundary voxels sit in the fine core.
grid = refine(coarse, sphere; factor=2, padding=1 // N)

geometry = Geometry([sphere], 0.0, AdaptiveAntiAliasing(SuperResolutionAntiAliasing(4)))
field = rasterize(geometry, grid)

# A dense array at the finest scale present in the grid, for plotting.
uniform = regrid(field)

# A field of region indices instead of rasterized values, regridded the same
# way, shows the tiling itself: each block of uniform color is one region.
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
