# VoxelShapes.jl

Place geometric shapes into a 3D voxel grid. Each shape carries a fill value or a gradient, and boundary voxels can be blended when a surface doesn't line up with the grid. `Array(world)` gives you a plain Julia array.

The library has no opinions about what the values mean. It works the same whether you're filling with electric susceptibilities, temperature fields, or binary masks.

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

world = World(
    (N, N, N),         # voxel grid dimensions
    (1/N, 1/N, 1/N),   # physical size of one voxel
    [sphere],
    0.0,               # background value
    NoAntiAliasing()
)

arr = Array(world)     # 64×64×64 Float64 array
```

## The `World` type

[`World`](@ref) is the central object. It holds the grid geometry, a tuple of shapes, the background fill value, and the anti-aliasing strategy. It's immutable, so use [`add_shape`](@ref) to build one up incrementally:

```julia
world = add_shape(world, FillableCylinder((0.5, 0.5, 0.5), 0.1, 0.4, 0.5))
```

Shapes are evaluated in order. The first shape whose containment test passes claims the voxel. If a shape evaluates to the background value, it's treated as transparent and evaluation falls through to the next shape. So you can punch holes by filling with the background value. The flip side: you can't deliberately paint the background value on top of a lower layer.

The `origin_xyz` keyword argument (default `(0, 0, 0)`) sets the world-space position of the front-bottom-left corner of voxel `(1, 1, 1)`.

```julia
world = World(
    (N, N, N),
    (1/N, 1/N, 1/N),
    shapes,
    0.0,
    NoAntiAliasing();
    origin_xyz = (-0.5, -0.5, -0.5)
)
```

## Examples

The [`examples/`](https://github.com/pvirally/VoxelShapes.jl/tree/main/examples) directory has a self-contained script for each feature. Run from the `examples/` folder:

```bash
julia --project=. 01_hello_sphere.jl
```

| Script | What it shows |
|---|---|
| `01_hello_sphere.jl` | Minimal example: one sphere, one world |
| `02_basic_shapes.jl` | All nine built-in primitives |
| `03_fills.jl` | Constant, radial, and axial gradients |
| `04_anti_aliasing.jl` | All five AA strategies side-by-side |
| `05_csg.jl` | Union, intersection, difference, hollow sphere |
| `06_rotation.jl` | Euler angles, axis-angle, explicit matrix |
| `07_showcase.jl` | Combined scene with multiple shapes and AA |

See the [API reference](@ref "API reference") for full signatures.
