# GPU

Rasterization runs on GPU with a single call:

```julia
using CUDA
arr = rasterize(geometry, region, CuArray)   # returns a CuArray{T, 3}
```

The kernel (`fill_voxel!`) is written with
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl)
and dispatches to `CUDABackend()` automatically. `rasterize(geometry,
region)` runs the same kernel on CPU. `rasterize(geometry, grid, CuArray)`
does the same over a `CompositeGrid`, region by region. It returns a
`CompositeField` whose `data` lives on the GPU.

## The `isbits` constraint

Only the `Geometry` reaches the kernel, so it has to be `isbits`. The
`Region` and `CompositeGrid` don't: they hold exact rationals and stay on
the host. Only the converted origin and voxel size of each region cross
over. This means:

- The shapes tuple and all shape fields must be stack-allocated (no heap
  references).
- Fill functions must be `isbits`. The built-in `ConstantFill`,
  `RadialGradient`, and `AxialGradient` all are. Closures that capture
  mutable objects are not.
- The background value must be `isbits` (plain numbers are fine).

The convenience constructors create closures (`_ -> fill_val`) that are
`isbits` only when `fill_val` is. Gradient fills need the inner struct
constructor (see [Fill functions](fills.md)).

## Performance notes

The kernel launches one GPU thread per voxel. That's efficient on large
grids. On small grids the launch overhead dominates, and CPU wins below
roughly 32³.

Shape evaluation is a type-stable recursive fold over the shape tuple.
Adding shapes costs compile time (more specializations), not runtime
branching: the compiler unrolls the tuple loop.

## Example

```julia
using VoxelShapes, CUDA

N = 256
dx = 1.0f0 / N   # Float32 for better GPU throughput

sphere = FillableSphere((0.5f0, 0.5f0, 0.5f0), 0.3f0, 1.0f0)
geometry  = Geometry([sphere], 0.0f0, SubpixelAntiAliasing())
region = Region((N, N, N), (1 // N, 1 // N, 1 // N), (1 // 2, 1 // 2, 1 // 2))

arr = rasterize(geometry, region, CuArray; precision=Float32)   # runs on GPU
```

Use `Float32` where you can. It packs twice as densely as `Float64` in GPU
SIMD, and it's usually accurate enough for voxel geometry.
