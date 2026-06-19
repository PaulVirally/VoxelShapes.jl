# 01_hello_sphere.jl
#
# The minimal VoxelShapes example. One sphere, one World, one array.
#
# A World describes the voxel grid and the shapes in it. Calling Array(world)
# runs the rasterizer and returns a plain Julia array you can slice and plot.
#
# Run from this directory:
#   julia --project=. 01_hello_sphere.jl

using VoxelShapes
using GLMakie

N = 64
dx = 1.0 / N

sphere = FillableSphere((0.5, 0.5, 0.5), 0.3, 1.0)

world = World(
    (N, N, N),
    (dx, dx, dx),
    [sphere],
    0.0,            # background value
    NoAntiAliasing()
)

arr = Array(world)

fig = Figure(size = (500, 500))
ax  = Axis(fig[1, 1], title = "Sphere, center z-slice", aspect = DataAspect())
heatmap!(ax, arr[:, :, N ÷ 2], colormap = :grays, colorrange = (0, 1))
hidedecorations!(ax, ticks = false, ticklabels = false)
save(joinpath(@__DIR__, "01_hello_sphere.png"), fig)
display(fig)
