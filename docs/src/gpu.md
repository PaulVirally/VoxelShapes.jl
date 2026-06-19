# GPU

Rasterization runs on GPU with a single call:

```julia
using CUDA
arr = CuArray(world)   # returns a CuArray{T, 3}
```

The kernel (`fill_voxel!`) is written with [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) and dispatches to `CUDABackend()` automatically. `Array(world)` runs the same kernel on CPU.

## The `isbits` constraint

For a `World` to be GPU-compatible, every field must be `isbits`. This means:

- The shapes tuple and all shape fields must be stack-allocated (no heap references).
- Fill functions must be `isbits`. The built-in `ConstantFill`, `RadialGradient`, and `AxialGradient` all are. Closures that capture mutable objects are not.
- The background value and voxel size must be `isbits` (plain numbers are fine).

The convenience constructors create closures (`_ -> fill_val`) that are `isbits` only for `isbits` `fill_val`. Gradient fills require the inner struct constructor (see [Fill functions](fills.md)).

## Performance notes

The kernel launches one GPU thread per voxel. That's efficient on large grids. On small grids the launch overhead dominates, and CPU wins below roughly 32³.

Shape evaluation is a type-stable recursive fold over the shape tuple. Adding shapes costs compile time (more specializations) but not runtime branching, since the compiler unrolls the tuple loop.

## Example

```julia
using VoxelShapes, CUDA

N = 256
dx = 1.0f0 / N   # Float32 for better GPU throughput

sphere = FillableSphere((0.5f0, 0.5f0, 0.5f0), 0.3f0, 1.0f0)
world  = World((N, N, N), (dx, dx, dx), [sphere], 0.0f0, SubpixelAntiAliasing())

arr = CuArray(world)   # runs on GPU
```

Use `Float32` where you can. It packs twice as densely as `Float64` in GPU SIMD and is usually accurate enough for voxel geometry.
