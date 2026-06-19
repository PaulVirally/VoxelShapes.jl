# 06_rotation.jl
#
# The Rotated wrapper tilts any shape in world space. The underlying shape
# stays axis-aligned; Rotated maps each query point into local frame before
# the containment test. This means every shape gets rotation for free.
#
# Three ways to specify the rotation:
#   Rotated(shape, (αx, αy, αz))        intrinsic ZYX Euler angles (radians)
#   Rotated(shape, axis, angle)          axis-angle
#   Rotated(shape, R)                    explicit 3×3 SMatrix
#
# The pivot defaults to center(shape). You can pass an explicit pivot as
# Rotated(shape, R, pivot) or Rotated(shape, angles, pivot).

using VoxelShapes
using GLMakie

N  = 80
dx = 1.0 / N
c  = (0.5, 0.5, 0.5)
aa = SubpixelAntiAliasing()

# A rectangular cuboid, easy to see when rotated.
box = FillableCuboid(c, (0.55, 0.25, 0.2), 1.0)

rotations = [
    ("No rotation",         box),
    ("45° around Z\n(Euler angles)",  Rotated(box, (0.0, 0.0, π/4))),
    ("30° around Y\n(axis-angle)",    Rotated(box, (0.0, 1.0, 0.0), π/6)),
    ("45° around Z,\n30° around X",  Rotated(box, (π/6, 0.0, π/4))),
]

fig = Figure(size = (900, 260))
for (i, (label, shape)) in enumerate(rotations)
    world = World((N, N, N), (dx, dx, dx), [shape], 0.0, aa)
    arr   = Array(world)
    ax    = Axis(fig[1, i], title = label, aspect = DataAspect())
    heatmap!(ax, arr[:, :, N ÷ 2], colormap = :grays, colorrange = (0, 1))
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "06_rotation.png"), fig)
display(fig)
