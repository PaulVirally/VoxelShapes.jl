# 03_fills.jl
#
# A fill function sets the value a voxel gets once it's inside a shape:
# ConstantFill, RadialGradient, or AxialGradient. Gradients need the inner
# struct constructor; the convenience constructors only take a constant
# value.

using VoxelShapes
using StaticArrays: SVector
using GLMakie

N = 64
region = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))
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
# Axis 3 is the z-axis for an ellipsoid. Local z goes from -1 (bottom) to +1
# (top).
f_axial = AxialGradient(3, 0.0, 1.0)  # dark bottom, bright top
sphere_axial = FillableEllipsoid{Float64,typeof(f_axial),typeof(interp)}(c, r, f_axial, interp)

# AxialGradient on a cylinder: axis 2 is the axial fraction along the tube.
cyl_center = SVector{3,Float64}(0.5, 0.5, 0.5)
f_cyl = AxialGradient(2, 0.0, 1.0)
cylinder_axial = FillableCylinder{Float64,typeof(f_cyl),typeof(interp)}(
    cyl_center, 0.25, 0.35, 3, f_cyl, interp
)

aa = SuperResolutionAntiAliasing(4)
geometry_for(shape) = Geometry([shape], 0.0, aa)

panels = [
    ("Constant (1.0)",       geometry_for(sphere_const),    arr -> arr[:, :, N ÷ 2]),
    ("Radial gradient",      geometry_for(sphere_radial),   arr -> arr[:, :, N ÷ 2]),
    ("Axial gradient (sphere)", geometry_for(sphere_axial), arr -> arr[:, :, N ÷ 2]),
    ("Axial gradient (cylinder)", geometry_for(cylinder_axial), arr -> arr[:, :, N ÷ 2]),
]

fig = Figure(size = (900, 280))
for (i, (title, geometry, slice_fn)) in enumerate(panels)
    arr = rasterize(geometry, region)
    ax  = Axis(fig[1, i], title = title, aspect = DataAspect())
    heatmap!(ax, slice_fn(arr), colormap = :inferno, colorrange = (0, 1))
    hidedecorations!(ax)
end
save(joinpath(@__DIR__, "03_fills.png"), fig)
display(fig)
