# 02_basic_shapes.jl
#
# A tour of every built-in shape. Each one is rasterized on its own so they
# can sit side-by-side. The z = 0.5 slice works for all of them except the
# capsule, which points along z, so that one gets a y-slice instead.
#
# Shapes covered: Sphere, Ellipsoid, Cuboid, Cylinder, Torus, Capsule,
# Cone (frustum), HalfSpace, Slab.

using VoxelShapes
using GLMakie

N = 64
aa = NoAntiAliasing()
region = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))

geometry_for(shape) = Geometry([shape], 0.0, aa)

c = (0.5, 0.5, 0.5)  # shared center

shapes = [
    ("Sphere",     FillableSphere(c, 0.3, 1.0)),
    ("Ellipsoid",  FillableEllipsoid(c, (0.35, 0.2, 0.25), 1.0)),
    ("Cuboid",     FillableCuboid(c, (0.6, 0.4, 0.5), 1.0)),
    ("Cylinder",   FillableCylinder(c, 0.25, 0.35, 1.0)),
    ("Torus",      FillableTorus(c, 0.28, 0.09, 1.0)),
    ("Capsule",    FillableCapsule((0.25, 0.5, 0.5), (0.75, 0.5, 0.5), 0.15, 1.0)),
    ("Cone",       FillableCone(c, 0.32, 0.0, 0.35, 1.0)),
    ("HalfSpace",  FillableHalfSpace(c, (1.0, 0.0, 0.0), 1.0)),
    ("Slab",       FillableSlab(c, (0.0, 0.0, 1.0), 0.15, 1.0)),
]

fig = Figure(size = (960, 640))
for (i, (name, shape)) in enumerate(shapes)
    arr   = rasterize(geometry_for(shape), region)
    row   = (i - 1) ÷ 3 + 1
    col   = (i - 1) % 3 + 1
    slice = name == "Capsule" ? arr[:, N ÷ 2, :] : arr[:, :, N ÷ 2]
    ax    = Axis(fig[row, col], title = name, aspect = DataAspect())
    heatmap!(ax, slice, colormap = :grays, colorrange = (0, 1))
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "02_basic_shapes.png"), fig)
display(fig)
