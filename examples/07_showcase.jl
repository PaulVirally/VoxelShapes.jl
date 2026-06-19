# 07_showcase.jl
#
# An "atom on a stand" scene that exercises the full library surface:
#
#   Shapes:     FillableSphere, FillableEllipsoid, FillableCylinder,
#               FillableCuboid, FillableCone, FillableTorus, FillableSlab,
#               FillableCapsule
#   Fills:      RadialGradient (nucleus), AxialGradient (pedestal),
#               constant values (everything else)
#   CSG:        csg_diff for a hollow nucleus and a hollow electron
#   Transforms: Rotated for the polar orbital ring
#   AA:         AdaptiveAntiAliasing wrapping SuperResolutionAntiAliasing
#
# Scene layout (z = 0 bottom, z = 1 top):
#   z ≈ 0.00 to 0.10   ground slab
#   z ≈ 0.10 to 0.24   base cuboid
#   z ≈ 0.24 to 0.50   pedestal cylinder  (axial gradient: dim at base, bright at top)
#   z ≈ 0.50 to 0.56   antenna capsule
#   z ≈ 0.60        nucleus  (hollow sphere with radial gradient fill)
#   z ≈ 0.60        equatorial ring (torus, xy-plane) and polar ring (torus, xz-plane)
#   z ≈ 0.60        electrons (small spheres; one is hollow via CSG)
#   z ≈ 0.78 to 0.92   top cone antenna

using VoxelShapes
using StaticArrays: SVector
using GLMakie

N  = 128
dx = 1.0 / N

interp = LinearInterpolation()

# ── Gradient fills require the inner struct constructor ───────────────────────

# Nucleus shell: bright core (r=0), dim at the surface (r=1)
rg = RadialGradient(1.0, 0.25)
nucleus_shell = FillableEllipsoid{Float64,typeof(rg),typeof(interp)}(
    SVector{3,Float64}(0.5, 0.5, 0.60),
    SVector{3,Float64}(0.14, 0.14, 0.14),
    rg, interp
)
# Hollow it out with CSG
nucleus_cavity = FillableSphere((0.5, 0.5, 0.60), 0.10, 1.0)
nucleus = csg_diff(nucleus_shell, nucleus_cavity)

# Pedestal: axis 2 in local coords is the axial fraction ∈ [-1, +1]
ag = AxialGradient(2, 0.15, 0.75)
pedestal = FillableCylinder{Float64,typeof(ag),typeof(interp)}(
    SVector{3,Float64}(0.5, 0.5, 0.37),
    0.055,   # radius
    0.13,    # half-height  (spans z ≈ 0.24 to 0.50)
    3,       # z axis
    ag, interp
)

# ── Constant-fill shapes ──────────────────────────────────────────────────────

ground      = FillableSlab((0.5, 0.5, 0.05), (0.0, 0.0, 1.0), 0.05, 0.2)
base_box    = FillableCuboid((0.5, 0.5, 0.17), (0.38, 0.38, 0.07), 0.35)
antenna_cap = FillableCapsule((0.5, 0.5, 0.50), (0.5, 0.5, 0.56), 0.013, 0.6)

# Equatorial ring: torus in the xy-plane (axis = 3)
eq_ring = FillableTorus((0.5, 0.5, 0.60), 0.26, 0.033, 0.7)

# Polar ring: same torus rotated 90° around x → lies in the xz-plane
pol_ring = Rotated(
    FillableTorus((0.5, 0.5, 0.60), 0.26, 0.033, 0.5),
    (π/2, 0.0, 0.0)
)

# Electrons at cardinal points of both rings
e_r = 0.26
electrons_solid = [
    FillableSphere((0.5 + e_r, 0.5,      0.60),       0.033, 0.90),  # equatorial +x
    FillableSphere((0.5,       0.5 + e_r, 0.60),       0.033, 0.90),  # equatorial +y
    FillableSphere((0.5,       0.5,       0.60 + e_r), 0.033, 0.75),  # polar +z
    FillableSphere((0.5,       0.5,       0.60 - e_r), 0.033, 0.75),  # polar -z
]

# One hollow electron (CSG) on the -x equatorial position
e_outer = FillableSphere((0.5 - e_r, 0.5, 0.60), 0.033, 0.90)
e_inner = FillableSphere((0.5 - e_r, 0.5, 0.60), 0.018, 0.90)
hollow_electron = csg_diff(e_outer, e_inner)

# Top cone antenna
top_cone = FillableCone((0.5, 0.5, 0.85), 0.020, 0.0, 0.07, 0.95)

# ── Assemble world ────────────────────────────────────────────────────────────

all_shapes = [
    ground, base_box, pedestal, antenna_cap,
    nucleus,
    eq_ring, pol_ring,
    electrons_solid...,
    hollow_electron,
    top_cone,
]

aa    = AdaptiveAntiAliasing(SuperResolutionAntiAliasing(4))
world = World((N, N, N), (dx, dx, dx), all_shapes, 0.0, aa)
arr   = Array(world)

# ── Visualize ────────────────────────────────────────────────────────────────
# Left:  3D isosurface rendered with GLMakie volume()
# Right: three orthogonal slices through the scene

iz_nucleus = round(Int, 0.60 / dx)
iy_center  = round(Int, 0.50 / dx)
ix_center  = round(Int, 0.50 / dx)

slice_specs = [
    ("XY  z ≈ 0.60  (nucleus level)",     arr[:, :, iz_nucleus]'),
    ("XZ  y ≈ 0.50  (equatorial cut)",    arr[:, iy_center, :]'),
    ("YZ  x ≈ 0.50  (polar cut)",         arr[ix_center, :, :]'),
]

fig = Figure(size = (1300, 750))

ax3d = Axis3(fig[1:3, 1], title = "3D isosurface (isovalue = 0.15)",
             aspect = :data, azimuth = 1.25π, elevation = 0.15π)
volume!(ax3d, 0.0..1.0, 0.0..1.0, 0.0..1.0, arr,
        algorithm = :iso, isovalue = 0.15,
        colormap = :linear_wcmr_100_45_c42_n256, transparency = false)

for (row, (title, slice)) in enumerate(slice_specs)
    ax = Axis(fig[row, 2], title = title, aspect = DataAspect())
    heatmap!(ax, slice, colormap = :linear_wcmr_100_45_c42_n256, colorrange = (0, 1))
    hidedecorations!(ax)
end

save(joinpath(@__DIR__, "07_showcase.png"), fig)
display(fig)
