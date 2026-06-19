# 05_csg.jl
#
# Constructive Solid Geometry (CSG) lets you combine shapes using boolean
# operations:
#
#   csg_union(a, b)      a OR b  (voxel in either shape)
#   csg_intersect(a, b)  a AND b (voxel in both)
#   csg_diff(a, b)       a AND NOT b (a with b carved out)
#   csg_complement(a)    NOT a (everything outside a; rarely rendered alone)
#
# Fill always delegates to the first (left) operand. You can chain CSG
# operations to build complex shapes from simple primitives.

using VoxelShapes
using GLMakie

N = 80
dx = 1.0 / N
aa = SubpixelAntiAliasing()

sphere = FillableSphere((0.45, 0.5, 0.5), 0.28, 1.0)
cube   = FillableCuboid((0.55, 0.5, 0.5), (0.45, 0.45, 0.45), 0.8)

# Hollow sphere: full sphere minus a slightly smaller concentric sphere
outer_shell  = FillableSphere((0.5, 0.5, 0.5), 0.3, 1.0)
inner_cavity = FillableSphere((0.5, 0.5, 0.5), 0.22, 1.0)

panels = [
    ("Union",        csg_union(sphere, cube)),
    ("Intersection", csg_intersect(sphere, cube)),
    ("Difference",   csg_diff(sphere, cube)),
    ("Hollow sphere\n(diff of two spheres)", csg_diff(outer_shell, inner_cavity)),
]

fig = Figure(size = (900, 260))
for (i, (label, shape)) in enumerate(panels)
    world = World((N, N, N), (dx, dx, dx), [shape], 0.0, aa)
    arr   = Array(world)
    ax    = Axis(fig[1, i], title = label, aspect = DataAspect())
    heatmap!(ax, arr[:, :, N ÷ 2], colormap = :thermal, colorrange = (0, 1))
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "05_csg.png"), fig)
display(fig)
