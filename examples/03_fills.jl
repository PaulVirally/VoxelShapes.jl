# 03_fills.jl
#
# Fill functions control what value a voxel gets once it's inside a shape.
# The default is a constant (the fill_val argument to every shape constructor).
# ConstantFill, RadialGradient, and AxialGradient let you vary the value
# spatially across the shape's volume.
#
# Because the convenience constructors only accept a constant value, gradient
# fills require passing the fill function directly via the inner struct
# constructor. StaticArrays.SVector is used for the center and radii fields.

using VoxelShapes
using StaticArrays: SVector
using GLMakie

N = 64
dx = 1.0 / N
c  = SVector{3,Float64}(0.5, 0.5, 0.5)
r  = SVector{3,Float64}(0.3, 0.3, 0.3)

interp = LinearInterpolation()

# ConstantFill: every voxel inside gets the same value.
f_const = ConstantFill(1.0)
sphere_const = FillableEllipsoid{Float64,typeof(f_const),typeof(interp)}(c, r, f_const, interp)

# RadialGradient: interpolates from inner_value at the center (r=0)
# to outer_value at the surface (r=1).
f_radial = RadialGradient(1.0, 0.0)   # bright core, dark shell
sphere_radial = FillableEllipsoid{Float64,typeof(f_radial),typeof(interp)}(c, r, f_radial, interp)

# AxialGradient: interpolates along one local axis of the shape.
# Axis 3 is the z-axis for an ellipsoid; local z goes from -1 (bottom) to +1 (top).
f_axial = AxialGradient(3, 0.0, 1.0)  # dark bottom, bright top
sphere_axial = FillableEllipsoid{Float64,typeof(f_axial),typeof(interp)}(c, r, f_axial, interp)

# AxialGradient on a cylinder: axis 2 is the axial fraction along the tube.
cyl_center = SVector{3,Float64}(0.5, 0.5, 0.5)
f_cyl = AxialGradient(2, 0.0, 1.0)
cylinder_axial = FillableCylinder{Float64,typeof(f_cyl),typeof(interp)}(
    cyl_center, 0.25, 0.35, 3, f_cyl, interp
)

aa = SuperResolutionAntiAliasing(4)
world_for(shape) = World((N, N, N), (dx, dx, dx), [shape], 0.0, aa)

panels = [
    ("Constant (1.0)",       world_for(sphere_const),    arr -> arr[:, :, N ÷ 2]),
    ("Radial gradient",      world_for(sphere_radial),   arr -> arr[:, :, N ÷ 2]),
    ("Axial gradient (sphere)", world_for(sphere_axial), arr -> arr[:, :, N ÷ 2]),
    ("Axial gradient (cylinder)", world_for(cylinder_axial), arr -> arr[:, :, N ÷ 2]),
]

fig = Figure(size = (900, 280))
for (i, (title, world, slice_fn)) in enumerate(panels)
    arr = Array(world)
    ax  = Axis(fig[1, i], title = title, aspect = DataAspect())
    heatmap!(ax, slice_fn(arr), colormap = :inferno, colorrange = (0, 1))
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "03_fills.png"), fig)
display(fig)
