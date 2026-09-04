# VoxelShapes.jl

[![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://paulvirally.github.io/VoxelShapes.jl/stable)
[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://paulvirally.github.io/VoxelShapes.jl/dev)
[![CI](https://github.com/PaulVirally/VoxelShapes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PaulVirally/VoxelShapes.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/PaulVirally/VoxelShapes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PaulVirally/VoxelShapes.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Version](https://juliahub.com/docs/General/VoxelShapes/stable/version.svg)](https://juliahub.com/ui/Packages/General/VoxelShapes)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

![A geometry built from every available shape primitive](examples/07_showcase.png)

Place geometric shapes into a 3D voxel grid. Each shape carries a fill value (or a gradient function), and you can pick how boundary voxels are blended when a surface doesn't line up with the grid. Call `rasterize(geometry, region)` and you get a plain Julia array.

The main use case is any grid-based simulation that needs per-voxel scalar or vector properties defined by geometry. The library has no opinions about what the values mean.

## Installation

```julia
using Pkg
Pkg.add("VoxelShapes")
```

## Quick start

```julia
using VoxelShapes

N = 64
sphere = FillableSphere((0.5, 0.5, 0.5), 0.3, 1.0)

geometry  = Geometry([sphere], 0.0, NoAntiAliasing())    # what to draw
region = Region((N, N, N), (1//N, 1//N, 1//N), (1//2, 1//2, 1//2))  # where to sample

arr = rasterize(geometry, region)    # returns a 64×64×64 Float64 array
```

A `Geometry` holds the shapes, the background value, and the anti-aliasing strategy; it says nothing about resolution, so the same geometry can be rasterized at any scale. A `Region` holds the cell count, voxel size, and center of a single block of uniform voxels, in exact rational arithmetic. Both are immutable. To add more shapes, use `add_shape(geometry, shape)`, which returns a new `Geometry` with the shape appended.

![Center z-slice of a sphere rasterized at 64³](examples/01_hello_sphere.png)

## Shapes

Nine primitives are included.

![Center slices of all nine built-in shape primitives](examples/02_basic_shapes.png)

```julia
FillableSphere(center, radius, fill_val)
FillableEllipsoid(center, (rx, ry, rz), fill_val)
FillableCuboid(center, (lx, ly, lz), fill_val)      # lx/ly/lz are full side lengths
FillableCylinder(center, radius, half_height, fill_val; axis=3)
FillableTorus(center, major_radius, minor_radius, fill_val; axis=3)
FillableCapsule(point_a, point_b, radius, fill_val)
FillableCone(center, base_radius, top_radius, half_height, fill_val; axis=3)
FillableSlab(point, normal, half_thickness, fill_val)
FillableHalfSpace(point, normal, fill_val)
```

`FillableCone` with `top_radius = 0` is a true cone; unequal nonzero radii give a frustum.

Shapes are evaluated in the order they were added to the geometry. The first shape whose containment test passes claims the voxel, with one exception: the background value doubles as a transparency sentinel. A shape that covers a voxel but produces a value equal to the background is treated as transparent there, so the next shape (or the background) shows through. This lets you punch holes by filling with the background value, but it also means you cannot deliberately paint the background value on top of a lower layer.

## Fill functions

The fill value argument can be any constant, but if you want spatial variation you can pass one of the built-in fill structs instead.

![Constant fill, radial gradient, and axial gradient on a sphere and cylinder](examples/03_fills.png)

```julia
# Radial gradient: interpolates from inner_value at the center to outer_value at the surface
f = RadialGradient(1.0, 0.0)   # bright core, transparent shell

# Axial gradient: interpolates along a local axis from -1 to +1
f = AxialGradient(3, 0.0, 1.0) # dark bottom, bright top (local z-axis)
```

Because the convenience constructors wrap their `fill_val` argument in a closure, using a gradient requires the inner struct constructor. See `examples/03_fills.jl` for the full syntax.

Any callable that takes a 3-tuple of local coordinates and returns a scalar works as a fill function, as long as it is `isbits`-compatible (required for GPU use).

## CSG

Shapes can be combined with boolean operations.

![Union, intersection, difference, and a hollow sphere](examples/05_csg.png)

```julia
csg_union(a, b)       # inside a or b
csg_intersect(a, b)   # inside both a and b
csg_diff(a, b)        # inside a but not b
csg_complement(a)     # everything outside a
```

Fill always delegates to the first operand. Operations can be nested.

## Rotation

`Rotated` wraps any shape and maps query points into local frame before the containment test. This means every shape gets rotation without needing its own rotation logic.

![A rectangular box at four rotation angles](examples/06_rotation.png)

```julia
Rotated(shape, (αx, αy, αz))     # intrinsic ZYX Euler angles in radians
Rotated(shape, axis, angle)       # axis-angle
Rotated(shape, R)                 # explicit 3×3 SMatrix
```

The pivot defaults to `center(shape)`. You can pass an explicit pivot as the last argument.

## Anti-aliasing

When a surface doesn't align with the voxel grid, boundary voxels need some treatment.

![Five anti-aliasing strategies compared at full resolution and zoomed in on a sphere edge](examples/04_anti_aliasing.png)

`NoAntiAliasing` is a hard point test at the voxel center. `SuperResolutionAntiAliasing(n)` divides each voxel into n³ sub-samples and averages them. `SubpixelAntiAliasing` uses the signed distance function to estimate coverage analytically, one evaluation per voxel, with no inner loop. `GaussianAntiAliasing(σ, kernel_size)` convolves the boundary with a Gaussian for a softer edge. `AdaptiveAntiAliasing(inner)` wraps any strategy and skips the stencil for voxels that are clearly inside or outside.

For most use cases `AdaptiveAntiAliasing(SuperResolutionAntiAliasing(4))` is a good starting point.

## Interpolation

When sub-voxel samples are averaged during anti-aliasing, the interpolation strategy controls how they are combined. Each shape carries its own, set at construction time via the `interpolation` keyword argument.

| Strategy | Description |
|---|---|
| `LinearInterpolation()` | Weighted arithmetic mean. The default. |
| `HarmonicInterpolation()` | Weighted harmonic mean. For positive-definite quantities. |
| `GeometricMeanInterpolation()` | Weighted geometric mean. Also for positive-definite quantities. |
| `MaxInterpolation()` | Maximum sample value, ignoring weights. |
| `MinInterpolation()` | Minimum sample value, ignoring weights. |
| `DielectricInterpolation()` | Linear interpolation of electric susceptibility χ. |
| `MetalInterpolation()` | Interpolates via the complex refractive index, then recovers χ. |

The last two implement interpolation schemes from computational electrodynamics. They are functionally identical to `LinearInterpolation` and `HarmonicInterpolation` except for how the intermediate values are computed.

## GPU

```julia
using CUDA
arr = rasterize(geometry, region, CuArray)   # runs the rasterization kernel on the GPU
```

The geometry, all shapes, and all fill functions must be `isbits`-compatible. Everything built into VoxelShapes is. Custom fill functions must also be `isbits`: no closures that capture heap-allocated objects.

## Grids and refinement

A `Region` is one block of uniform voxels. Real problems often want higher resolution near a surface and coarser voxels away from it, without paying for a uniform fine grid everywhere. A `CompositeGrid` glues several `Region`s of different resolution into one grid that still rasterizes as a single object.

```julia
grid = refine(Region((8, 8, 8), (1//16, 1//16, 1//16)),
              ((0//1, 0//1, 0//1), (1//8, 1//8, 1//8)))
```

This carves the region into a refined core plus the six leftover slabs that fill out the rest of it, seven regions in total, matching the reference case from GilaElectromagnetics.jl's own composite-volume tests. `refine(grid, shape; factor=2, padding=0)` does the same starting from a shape's `bounding_box` instead of a box you specify by hand, so refining around a sphere is one call:

```julia
sphere = FillableSphere((0.0, 0.0, 0.0), 0.05, 1.0)
grid   = refine(Region((8, 8, 8), (1//16, 1//16, 1//16)), sphere; factor=2, padding=1//16)
```

`rasterize(geometry, grid)` returns a `CompositeField`: one flat vector, one entry per voxel, region by region in the order `refine` produced them. `vec(field)` is that flat vector, and `collect(eachregion(field))` gives one 3D array per region; both line up with `regionview` and `eachregion`. `regrid(field)` resamples the whole thing onto one uniform array at the finest scale present, which is what you want for a quick `heatmap!`, at the cost of the memory a uniform grid at that resolution needs. See `examples/08_refine.jl` and the [Grids and refinement](https://paulvirally.github.io/VoxelShapes.jl/stable/grids/) docs page for the full picture, including a slice colored by region index.

This whole grid representation, and `refine` itself, mirrors [GilaElectromagnetics.jl](https://github.com/PaulVirally/GilaElectromagnetics.jl)'s own composite volumes, so a `CompositeField` rasterized here lines up with Gila's degrees of freedom without any reshaping.

## Extending

To implement a custom shape, define:

```julia
Base.in(point::NTuple{3,T}, shape::MyShape)                            # containment test
Base.fill(shape::MyShape, voxel_center::NTuple{3,T}, voxel_size::NTuple{3,T})  # fill value
VoxelShapes.interpolation(shape::MyShape)                              # interpolation strategy
VoxelShapes.sdf(shape::MyShape, point::NTuple{3,T})                    # signed distance function
```

`has_exact_sdf` defaults to `false`. Set it to `true` if your SDF is a true Euclidean distance, which lets `AdaptiveAntiAliasing` skip boundary checks for interior and exterior voxels.

Optionally, define `VoxelShapes.bounding_box(shape::MyShape)` to return a conservative axis-aligned `(lower, upper)` box (infinite extents as `±Inf`). The default is all-infinite, so `refine(grid, shape)` on a shape without this method refines the whole grid instead of just the neighborhood of the shape.
