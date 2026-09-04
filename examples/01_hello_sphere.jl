# 01_hello_sphere.jl
#
# The simplest VoxelShapes example: one sphere, one Geometry, one Region.
#
# A Geometry says what to draw: the shapes, a background value, and an
# anti-aliasing strategy. A Region says where to sample: a voxel grid.
# rasterize(geometry, region) runs the rasterizer and returns a plain Julia
# array you can slice and plot.
#
# Run from this directory:
#   julia --project=. 01_hello_sphere.jl

using VoxelShapes
using GLMakie

N = 64

sphere = FillableSphere((0.5, 0.5, 0.5), 0.3, 1.0)

geometry = Geometry(
    [sphere],
    0.0,            # background value
    NoAntiAliasing()
)
region = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))

arr = rasterize(geometry, region)

fig = Figure(size = (500, 500))
ax  = Axis(fig[1, 1], title = "Sphere, center z-slice", aspect = DataAspect())
heatmap!(ax, arr[:, :, N ÷ 2], colormap = :grays, colorrange = (0, 1))
hidedecorations!(ax, ticks = false, ticklabels = false)
save(joinpath(@__DIR__, "01_hello_sphere.png"), fig)
display(fig)
