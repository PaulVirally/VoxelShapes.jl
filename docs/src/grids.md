# Grids and refinement

A [`Region`](@ref) is a single block of uniform voxels: a cell count, a voxel
side length, and a center, all exact rationals. This page covers gluing
regions of different resolution into one [`CompositeGrid`](@ref), and using
`refine` to build one around a shape.

## Region

```julia
region = Region((64, 64, 64), (1//64, 1//64, 1//64), (1//2, 1//2, 1//2))
```

`region` spans `center - cells .* scale / 2` to `center + cells .* scale / 2`.
This one covers the unit cube. `cells`, `scale`, `center`, `lower_corner`,
`upper_corner`, and `voxel_centers` read the fields back out. Integer
arguments for `scale` and `center` are converted to `Rational{Int}`. `1` is a
valid cell count along any axis.

## Refining by a box

`refine(grid, box; factor=2, padding=0)` replaces every region the box
touches with a refined core plus the leftover slabs that fill out the rest of
that region.

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

Region 1 is the refined core, at half the parent scale. The rest are the
six complement slabs, in the fixed order xlo, xhi, ylo, yhi, zlo, zhi.
This order fixes the flat layout of any field over the grid. `refine` can
be nested or applied more than once. A `Region` passed directly to
`refine` is wrapped in a one-region `CompositeGrid` first.

## Refining by a shape

`refine` also takes a shape instead of a box:

```julia
grid = refine(coarse_region, sphere; factor=2, padding=1//64)
```

This refines around `bounding_box(sphere)`, grown by `padding` on every
side and clamped to the grid's own bounding box. Every built-in shape
implements `bounding_box`. A custom shape that doesn't gets the default
all-infinite box: refining around it refines the whole grid.

## Rasterizing into a composite grid

`rasterize(geometry, grid)` runs the same kernel `rasterize(geometry,
region)` does, once per region, and returns a [`CompositeField`](@ref): one
flat vector, one entry per voxel, in the region layout order above.

```julia
field = rasterize(geometry, grid)
vec(field)                    # flat Vector, length(grid) entries
collect(eachregion(field))    # one 3D array view per region, in layout order
```

`eachregion` and `regionview` are reshaped views into `vec(field)`, not
copies.

## `regrid` for plotting

A `CompositeField` isn't a dense array. There's nothing to hand `heatmap!`
directly. `regrid(field, scale=finest(field.grid))` resamples it onto one
uniform array covering the grid's bounding box, replicating each coarse
voxel over every finer output cell it contains:

```julia
uniform = regrid(field)          # at the finest scale present in the grid
heatmap(uniform[:, :, end÷2])
```

`regrid` runs on the CPU; a field rasterized on the GPU is copied to the host
first. See `examples/08_refine.jl` for a slice colored by region index.

See the [API reference](@ref "API reference") for full signatures.
