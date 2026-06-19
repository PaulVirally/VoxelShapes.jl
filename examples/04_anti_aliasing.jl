# 04_anti_aliasing.jl
#
# Anti-aliasing decides what happens to the voxels that straddle a shape's
# surface. Zoom in on a sphere edge to see the difference.
#
# NoAntiAliasing:              hard on/off at the voxel center.
# SuperResolutionAntiAliasing: subdivides each voxel into sub-samples,
#                              averages the results. Accurate but slow.
# SubpixelAntiAliasing:        uses the signed distance function (SDF) to
#                              estimate surface coverage analytically. O(1),
#                              GPU-safe.
# AdaptiveAntiAliasing:        skips the sub-sampling stencil for voxels
#                              clearly inside or outside; only boundary voxels
#                              pay the full cost. Wraps any other strategy.
# GaussianAntiAliasing:        convolves the boundary with a Gaussian kernel.
#                              Produces a soft edge rather than a crisp one.

using VoxelShapes
using GLMakie

N = 64
dx = 1.0 / N
sphere = FillableSphere((0.5, 0.5, 0.5), 0.3, 1.0)

strategies = [
    ("None",              NoAntiAliasing()),
    ("SuperRes(4)",       SuperResolutionAntiAliasing(4)),
    ("Subpixel",          SubpixelAntiAliasing()),
    ("Adaptive+SuperRes", AdaptiveAntiAliasing(SuperResolutionAntiAliasing(4))),
    ("Gaussian(σ=1, k=5)", GaussianAntiAliasing(1.0, 5)),
]

fig = Figure(size = (1100, 520))
for (i, (label, aa)) in enumerate(strategies)
    world = World((N, N, N), (dx, dx, dx), [sphere], 0.0, aa)
    arr   = Array(world)
    slice = arr[:, :, N ÷ 2]

    # Full slice
    ax1 = Axis(fig[1, i], title = label, aspect = DataAspect())
    heatmap!(ax1, slice, colormap = :grays, colorrange = (0, 1))
    hidedecorations!(ax1)

    # Zoomed edge (top-right quadrant of the sphere boundary)
    lo, hi = 3N ÷ 4 - 12, 3N ÷ 4 + 12
    ax2 = Axis(fig[2, i], aspect = DataAspect())
    heatmap!(ax2, slice[lo:hi, lo:hi], colormap = :grays, colorrange = (0, 1))
    hidedecorations!(ax2)
end

Label(fig[2, 1:end, Bottom()], "Bottom row: zoomed edge (top-right boundary region)")
save(joinpath(@__DIR__, "04_anti_aliasing.png"), fig)
display(fig)
