# VoxelShapes.jl

Place geometric shapes into a 3D voxel grid. Each shape carries a fill
value or a gradient. Boundary voxels can be blended when a surface doesn't
line up with the grid. `rasterize` runs the fill and gives you a plain
Julia array.

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

arr = rasterize(geometry, region)     # 64×64×64 Float64 array
```

## `Geometry` and `Region`

A [`Geometry`](@ref) says what to draw: a tuple of shapes, a background
value, and an anti-aliasing strategy. Resolution is set by the `Region`,
not the `Geometry`. It's immutable. Use [`add_shape`](@ref) to build one
up:

```julia
geometry = add_shape(geometry, FillableCylinder((0.5, 0.5, 0.5), 0.1, 0.4, 0.5))
```

Shapes are evaluated in order. The first shape whose containment test
passes claims the voxel. If a shape evaluates to the background value,
it's treated as transparent and evaluation falls through to the next
shape. You can punch holes by filling with the background value, but you
can't deliberately paint the background value on top of a lower layer.

A [`Region`](@ref) says where to sample: a cell count, a voxel side
length, and a center, all in exact rational arithmetic. It spans `center -
cells .* scale / 2` to `center + cells .* scale / 2`. A region centered on
`(1/2, 1/2, 1/2)` with `N` unit-fraction voxels per axis covers the unit
cube:

```julia
region = Region((N, N, N), (1//N, 1//N, 1//N), (1//2, 1//2, 1//2))
```

[`rasterize`](@ref) fills every voxel center of a region with the
geometry's value there. For regions of mixed resolution, see [Grids and
refinement](@ref).

## Examples

The [`examples/`](https://github.com/pvirally/VoxelShapes.jl/tree/main/examples)
directory has a self-contained script for each feature. Run from the
`examples/` folder:

```bash
julia --project=. 01_hello_sphere.jl
```

| Script | What it shows |
|---|---|
| `01_hello_sphere.jl` | Minimal example: one sphere, one geometry, one region |
| `02_basic_shapes.jl` | All nine built-in primitives |
| `03_fills.jl` | Constant, radial, and axial gradients |
| `04_anti_aliasing.jl` | All five AA strategies side-by-side |
| `05_csg.jl` | Union, intersection, difference, hollow sphere |
| `06_rotation.jl` | Euler angles, axis-angle, explicit matrix |
| `07_showcase.jl` | Combined geometry with multiple shapes and AA |
| `08_refine.jl` | A composite grid refined around a sphere, and `regrid` for plotting |

See the [API reference](@ref "API reference") for full signatures.
