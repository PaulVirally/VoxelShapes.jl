# Grids and refinement

A [`Region`](@ref) is a single block of uniform voxels: a cell count, a voxel side length, and a center, all exact rationals. Every earlier page rasterizes into one `Region`. This page covers gluing regions of different resolution together into one [`CompositeGrid`](@ref), and using `refine` to build one around a shape.

The exact-rational grid representation and `refine` come from GilaElectromagnetics.jl, the electromagnetic solver this library feeds geometry to. A `CompositeGrid` here reproduces the region layout of Gila's `GlaCmpVol` exactly, so a geometry rasterized into one lines up with Gila's degrees of freedom without any extra bookkeeping.

## Region

```julia
region = Region((64, 64, 64), (1//64, 1//64, 1//64), (1//2, 1//2, 1//2))
```

`region` spans `center - cells .* scale / 2` to `center + cells .* scale / 2`, so this one covers the unit cube. `cells`, `scale`, `center`, `lower_corner`, `upper_corner`, and `voxel_centers` read the fields back out. Integer arguments for `scale` and `center` are converted to `Rational{Int}`, and `1` is a valid cell count along any axis for a lone region.

## Refining by a box

`refine(grid, box; factor=2, padding=0)` replaces every region the box touches with a refined core plus the leftover slabs that fill out the rest of that region. Starting from a single 8x8x8 region at scale `1/16` and refining by 2 inside a box of side `1/8` centered on the origin gives the same seven-region layout Gila's reference test produces:

```julia
julia> grid = refine(Region((8, 8, 8), (1//16, 1//16, 1//16)),
                      ((0//1, 0//1, 0//1), (1//8, 1//8, 1//8)))
CompositeGrid (7 regions)
  [1] (8×8×8) cells, scale (1//32, 1//32, 1//32), center (0//1, 0//1, 0//1)
  [2] (2×8×8) cells, scale (1//16, 1//16, 1//16), center (-3//16, 0//1, 0//1)
  [3] (2×8×8) cells, scale (1//16, 1//16, 1//16), center (3//16, 0//1, 0//1)
  [4] (4×2×8) cells, scale (1//16, 1//16, 1//16), center (0//1, -3//16, 0//1)
  [5] (4×2×8) cells, scale (1//16, 1//16, 1//16), center (0//1, 3//16, 0//1)
  [6] (4×4×2) cells, scale (1//16, 1//16, 1//16), center (0//1, 0//1, -3//16)
  [7] (4×4×2) cells, scale (1//16, 1//16, 1//16), center (0//1, 0//1, 3//16)
```

Region 1 is the refined core, at half the parent scale; the rest are the six complement slabs, in the fixed order xlo, xhi, ylo, yhi, zlo, zhi. That order fixes the flat layout of any field over the grid: region 1's voxels come first, then region 2's, and so on. `refine` can be nested or applied more than once, each time carving further into whichever regions the next box touches.

A `Region` passed directly to `refine` is wrapped in a one-region `CompositeGrid` first.

## Refining by a shape

Building a box by hand is tedious, so `refine` also takes a shape:

```julia
grid = refine(coarse_region, sphere; factor=2, padding=1//64)
```

This refines around `bounding_box(sphere)`, grown by `padding` on every side and clamped to the grid's own bounding box. Every built-in shape implements `bounding_box`; a custom shape that doesn't gets the default all-infinite box, so refining around it just refines the whole grid. A padding of one coarse voxel is a reasonable default when the shape has an anti-aliased boundary, since the boundary voxels are the ones that benefit from the extra resolution.

## Rasterizing into a composite grid

`rasterize(geometry, grid)` runs the same kernel `rasterize(geometry, region)` does, once per region, and returns a [`CompositeField`](@ref): one flat vector, one entry per voxel, in the region layout order above.

```julia
field = rasterize(geometry, grid)
vec(field)                    # flat Vector, length(grid) entries
collect(eachregion(field))    # one 3D array view per region, in layout order
```

`vec(field)` is exactly the flat, per-cell susceptibility form Gila's composite operators accept, and `collect(eachregion(field))` is the per-region tensor form. Neither copies data: `eachregion` and `regionview` are reshaped views into `vec(field)`, so a Gila susceptibility array built from either sees whatever `rasterize` wrote.

## `regrid` for plotting

A `CompositeField` isn't a dense array, so there's nothing to hand `heatmap!` directly. `regrid(field, scale=finest(field.grid))` resamples it onto one uniform array covering the grid's bounding box, replicating each coarse voxel over every finer output cell it contains:

```julia
uniform = regrid(field)          # at the finest scale present in the grid
heatmap(uniform[:, :, end÷2])
```

This is the cheapest way to look at a composite field, at the cost of the memory a uniform array at that resolution needs. `regrid` runs on the CPU; a field rasterized on the GPU is copied to the host first. See `examples/08_refine.jl` for a complete script, including a slice colored by region index rather than by field value, which shows the tiling itself.

See the [API reference](@ref "API reference") for full signatures.
