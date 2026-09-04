# 01_hello_sphere.jl
#
# The minimal VoxelShapes example. One sphere, one Geometry, one Region.
#
# A Geometry describes what to draw: the shapes, the background value, and the
# anti-aliasing strategy. A Region describes where to sample: a voxel grid
# geometry. rasterize(geometry, region) runs the rasterizer and returns a plain
# Julia array you can slice and plot.
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
