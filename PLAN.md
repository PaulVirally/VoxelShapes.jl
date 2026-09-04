# PLAN: Composite grids and `refine` for VoxelShapes.jl (v0.2.0)

VoxelShapes.jl is the geometry front-end for GilaElectromagnetics.jl
(branch `paul-gila-operators`). Gila recently gained composite volumes
(`GlaCmpVol`) and a `refine` function that glues together regions of different
resolution. This plan brings the same concepts and the same interface to
VoxelShapes so that a geometry can be rasterized straight into a Gila composite
volume's layout.

This is a **breaking release (0.2.0)**. `World` is removed (not aliased) and
replaced by a `Geometry` (what to draw) plus a grid (`Region` or `CompositeGrid`,
where to sample).

## 1. Background: how Gila defines volumes

Reference source: `src/glaVol.jl`, `src/glaCmpVol.jl`, `src/glaFld.jl`,
`src/glaCmpOpr.jl` on the `paul-gila-operators` branch of
https://github.com/PaulVirally/GilaElectromagnetics.jl.

### `GlaVol`
```julia
struct GlaVol
    cel::NTuple{3,Integer}    # cells per axis
    scl::NTuple{3,Rational}   # cell side lengths, in wavelengths
    org::NTuple{3,Rational}   # CENTER of the volume, in wavelengths
    grd::Array{<:StepRange,1} # cell-center coordinates per axis
end
GlaVol(cel, scl, org=(0//1,0//1,0//1))
```
Cell centers: lower corner is `org - cel*scl/2`; center of cell `i` along an
axis is `lower + (i - 1/2)*scl`. Everything is exact `Rational` arithmetic.

Lower/upper corners used everywhere:
```julia
_lwrEdg(vol) = Tuple(first.(vol.grd) .- (vol.scl .// 2))   # == org - cel*scl//2
_uprEdg(vol) = Tuple(last.(vol.grd)  .+ (vol.scl .// 2))   # == org + cel*scl//2
```

### `GlaCmpVol`
`struct GlaCmpVol; regions::Vector{GlaVol}; end`. Invariants checked by the
constructor (`chkCmpVol`), each throwing `ArgumentError` with a descriptive
message:
1. Non-empty region list; every region has `cel >= 1` in every axis.
2. Every region has an **even** cell count in every axis.
3. Pairwise commensurate scales: `max(sclA,sclB) // min(sclA,sclB)` is an
   integer in every axis.
4. Disjoint interiors: `all(max.(lwrA,lwrB) .< min.(uprA,uprB))` must be false.
   Touching along faces/edges/corners is allowed.
5. Common grid: `(lwrA - lwrB) // gcd.(all scales)` is integer in every axis.
6. Partition parity (`chkParCmpVol`): with `maxScl = lcm.(sclA, sclB)`,
   `divA = maxScl ÷ sclA`, `divB = maxScl ÷ sclB`, require
   `mod.(celA, divA) == 0`, `mod.(celB, divB) == 0`, and
   `iseven.(celA ÷ divA .+ celB ÷ divB)` in every axis.
7. Exact tiling (`chkTilCmpVol`): sum of region volumes equals volume of
   bounding box of the union (exact Rational comparison).

Accessors: `regions(cvol)`, `nregions(cvol)`,
`finest(cvol)` (elementwise min scale), `coordinates(cvol)` (iterator yielding
`(pos::NTuple{3,Float64}, cellvolume::Float64, regionindex::Int)` per cell,
region by region, column-major within a region), `cellvolumes(cvol)`.
Also `Base.:(==)` and `Base.show`.

### `refine(cvol, box; factor=2, snap=:outward)`
- `box` is a `GlaVol` (its outer bounds are used) or a tuple
  `(center::NTuple{3,Rational}, dims::NTuple{3,Rational})` of center and full
  side lengths.
- `factor` is an `Integer` or `NTuple{3,Integer}`, every entry `>= 1`.
- `snap` must be `:outward` (only mode).
- For each region whose interior meets the box (`all(max.(boxLwr,regLwr) .< min.(boxUpr,regUpr))`):
  1. Intersect the box with the region.
  2. `_snpBox`: convert to zero-based, half-open cell index bounds:
     `rawLwr = floor((boxLwr - regLwr)/scl)`, `rawUpr = ceil((boxUpr - regLwr)/scl)`;
     then snap to even indices: `idxLwr = 2*fld(max(rawLwr,0),2)`,
     `idxUpr = 2*cld(min(rawUpr,cel),2)`.
  3. `_crvVol`: replace the region by the refined core
     (cells `(idxUpr-idxLwr)*fac`, scale `scl//fac`, same physical extent)
     followed by the non-empty complement slabs in the fixed order
     **xlo, xhi, ylo, yhi, zlo, zhi**:
     - xlo: cells `(0,0,0)`→`(idxLwr[1], cel[2], cel[3])`
     - xhi: `(idxUpr[1],0,0)`→`(cel[1], cel[2], cel[3])`
     - ylo: `(idxLwr[1],0,0)`→`(idxUpr[1], idxLwr[2], cel[3])`
     - yhi: `(idxLwr[1],idxUpr[2],0)`→`(idxUpr[1], cel[2], cel[3])`
     - zlo: `(idxLwr[1],idxLwr[2],0)`→`(idxUpr[1], idxUpr[2], idxLwr[3])`
     - zhi: `(idxLwr[1],idxLwr[2],idxUpr[3])`→`(idxUpr[1], idxUpr[2], cel[3])`
     Slabs keep the parent scale. A slab is dropped when it has zero cells.
- Regions the box misses are kept in place. The result is validated by the
  `GlaCmpVol` constructor. Refinement can be applied repeatedly and nested.
- `refine(vol::GlaVol, box; ...)` wraps into a one-region composite first.

Gila's reference test (`test/cmpVolTest.jl`): refining
`GlaVol((8,8,8), (1//16,1//16,1//16), (0,0,0))` with box
`((0,0,0), (1//8,1//8,1//8))` gives 7 regions: a core of `(8,8,8)` cells at
`1//32` centered at the origin, then slabs with cell counts
`(2,8,8), (2,8,8), (4,2,8), (4,2,8), (4,4,2), (4,4,2)` at `1//16`.
Region ordering defines the **flat degree-of-freedom layout** Gila uses for
fields and susceptibility, so VoxelShapes must reproduce it exactly.

### What Gila accepts as susceptibility over a composite volume (`_cmpSus`)
A `Number`; a function `pos::NTuple{3,Float64} -> scalar`; a
`Vector` of one 3-tensor per region (each sized `reg.cel`); a flat `Vector` with
one value per cell in region layout order (column-major within a region,
regions in order); or a flat vector with `3 * ncells` entries (DOF layout).
Over a single `GlaVol`: a 3-tensor sized `cel` or its `vec`. Susceptibility is
`ComplexF64`.

## 2. Target VoxelShapes API

### 2.1 Grids (`src/Grids.jl`, `module Grids`)
Exact Rational arithmetic throughout. All lengths are in the same units Gila
uses (wavelengths), but VoxelShapes is unit-agnostic and only documents this.

```julia
struct Region
    cells::NTuple{3,Int}
    scale::NTuple{3,Rational{Int}}    # voxel side lengths
    center::NTuple{3,Rational{Int}}   # center of the region
end
Region(cells::NTuple{3,Integer}, scale::NTuple{3,Rational}, center::NTuple{3,Rational}=(0//1,0//1,0//1))
# Also accept Integer entries in scale/center (convert to Rational{Int}).
# Validate cells >= 1 (ArgumentError). Even cell counts are NOT required for a
# lone Region (only for CompositeGrid, mirroring Gila).
```
Functions on `Region`:
- `cells(r)`, `scale(r)`, `Base.size(r) = cells(r)`, `Base.length(r) = prod(cells)`
- `VoxelShapes.center(r)` (extend `Types.center`), `lower_corner(r)`, `upper_corner(r)`
- `voxel_centers(r) -> NTuple{3, StepRange{Rational{Int},Rational{Int}}}`
  (same values as Gila `grd`)
- `Base.:(==)`, `Base.show`

```julia
struct CompositeGrid
    regions::Vector{Region}
    # inner constructor runs the seven invariant checks above; copies the vector
end
CompositeGrid(r::Region)
```
Functions on `CompositeGrid` (names identical to Gila):
- `regions`, `nregions`, `finest`, `coordinates` (iterator of
  `(pos::NTuple{3,Float64}, cellvolume::Float64, regionindex::Int)`, with
  `length`, `eltype`, `IteratorSize`), `cellvolumes` (**one entry per cell**, not
  three; VoxelShapes fields are scalar. Document this difference from Gila.)
- `bounding_box(g) -> (lower::NTuple{3,Rational{Int}}, upper::NTuple{3,Rational{Int}})`
- `Base.length(g)` = total cell count, `Base.:(==)`, `Base.show`
- `refine(g::CompositeGrid, box; factor=2, snap=:outward)` and
  `refine(r::Region, box; ...)` exactly as in Gila (section 1), where `box` is a
  `Region` or `(center, dims)` tuple of `NTuple{3,Rational}`.
- Shape overload (depends on section 2.2):
  `refine(g, shape::AbstractFillableShape; factor=2, padding=0, snap=:outward)`
  uses `bounding_box(shape)` widened by `padding` on every side, clamped to
  `bounding_box(g)`, and then refines. Bounds may be `Float64` (and `±Inf`);
  `_snpBox`'s floor/ceil handles floats after the clamp. If the clamped box has
  empty interior, return the grid unchanged.

Error messages should be as descriptive as Gila's (which region, which axis,
what numbers, how to fix).

### 2.2 Shape bounding boxes (`Types.jl` + every shape module)
New interface function in `Types`:
```julia
bounding_box(shape) -> (lower::NTuple{3,T}, upper::NTuple{3,T})
```
Axis-aligned, conservative (may be larger than the shape, never smaller).
Infinite extents use `±Inf`. Methods:
- Ellipsoid/Sphere: `center ± radii`
- Cuboid: `center ± half_lengths`
- Cylinder: `± half_height` on `axis`, `± radius` on the other two
- Cone: `± half_height` on `axis`, `± max(base_radius, top_radius)` otherwise
- Torus: `± minor_radius` on `axis`, `± (major_radius + minor_radius)` otherwise
- Capsule: `min.(a,b) - radius`, `max.(a,b) + radius`
- Slab: on axis `i`, if `abs(normal[i]) ≈ 1` then `point[i] ± half_thickness`,
  else `(-Inf, Inf)`
- HalfSpace: if `abs(normal[i]) ≈ 1`, bounded on one side (`normal[i] > 0` →
  upper bound `point[i]`, lower `-Inf`; `normal[i] < 0` the reverse); otherwise
  `(-Inf, Inf)`. Check the existing containment convention in `Slabs.jl` for
  which side is inside.
- Rotated: if any inner bound is infinite, return all-infinite. Otherwise map
  the 8 corners of the inner box to world space with the inverse of the
  world-to-local transform: `world = R' * (c - pivot) + pivot`, and take
  componentwise min/max.
- CSG: `UnionShape` = union of boxes; `IntersectionShape` = intersection
  (componentwise max of lowers, min of uppers); `DifferenceShape` = box of `a`;
  `ComplementShape` = all-infinite.
Add `bounding_box` to the "Extending" docs as an optional method (default:
all-infinite for unknown shapes, so refine-by-shape degrades to "refine
everything").

### 2.3 Geometry (`src/Geometry.jl`, included at top level, not a submodule)
```julia
struct Geometry{U, S<:Tuple, A<:AbstractAntiAliasing}
    shapes::S
    background::U
    aa::A
end
Geometry(shapes::Union{Tuple,AbstractVector}, background, aa::AbstractAntiAliasing=NoAntiAliasing())
add_shape(geometry, shape) -> Geometry
Base.eltype(geometry) = typeof(background)
sample(geometry, pos::NTuple{3,T}, voxel_size::NTuple{3,T}) -> U
```
`sample` is the existing `_fill_shapes` fold. It lets a geometry be handed to
Gila's function-of-position susceptibility form via
`pos -> sample(geometry, pos, vs)`.

### 2.4 Rasterization (`src/Rasterize.jl`, included at top level)
Single region:
```julia
rasterize(geometry, region::Region; precision::Type{<:AbstractFloat}=Float64) -> Array{U,3}
rasterize(geometry, region::Region, ::Type{CuArray}; precision=Float64) -> CuArray{U,3}
rasterize!(arr::Array{U,3}, geometry, region; precision=Float64)     # CPU backend
rasterize!(arr::CuArray{U,3}, geometry, region; precision=Float64)   # CUDABackend
```
The kernel is the existing `fill_voxel!` reworked to take the geometry plus a
lower corner and voxel size converted to `precision`
(`origin = precision.(lower_corner(region))`, `voxel_size = precision.(scale(region))`).
Voxel center of index `idx` is `origin + (idx - 1/2) * voxel_size`, which
matches Gila's cell centers. `precision` must match the shapes' `T`
(Float32 shapes need `precision=Float32`; the existing CUDA tests do this).

Composite:
```julia
struct CompositeField{U, V<:AbstractVector{U}}
    data::V                 # flat, one entry per cell, region layout order
    grid::CompositeGrid
    offsets::Vector{Int}    # offsets[i]+1 : offsets[i+1] is region i's block; length nregions+1
end
rasterize(geometry, grid::CompositeGrid; precision=Float64) -> CompositeField{U, Vector{U}}
rasterize(geometry, grid::CompositeGrid, ::Type{CuArray}; precision=Float64) -> CompositeField{U, CuVector{U}}
```
Implementation: allocate the flat buffer, and for each region run the
single-region kernel into `reshape(view(data, block), cells)` (one kernel launch
per region; KernelAbstractions accepts reshaped views on both backends. If a
backend rejects the view, rasterize into a temporary and `copyto!`).

Methods on `CompositeField`: `Base.vec(f) = f.data`, `Base.eltype`,
`Base.length`, `regions(f)`, `nregions(f)`, `regionview(f, i)` (reshape of the
view to `cells(region)`), `eachregion(f)`,
`regrid(f, scale=finest(f.grid)) -> Array{U,3}` (uniform array over the
bounding box; every uniform cell takes the value of the source cell containing
its center; errors like Gila if `scale` does not divide every region scale or
the bounding box; CPU only, copies GPU data to host first), `Base.show`,
`Base.:(==)`.

`vec(field)` is directly one of the forms Gila's `_cmpSus` accepts (per-cell
flat vector in layout order), and `collect(eachregion(field))` is the
per-region-tensor form.

### 2.5 Removal of `World`
Delete `World`, its constructors, `Base.Array(::World)`, `CUDA.CuArray(::World)`,
`Base.fill!(::Array, ::World)`, `Base.similar(::World)`, `idx2global_xyz`.
Do not keep an alias. Update every reference (see section 4).

### 2.6 Exports (top-level `VoxelShapes` module)
Add: `Geometry, sample, rasterize, rasterize!, CompositeField, regionview,
eachregion, regrid, Region, CompositeGrid, refine, regions, nregions, finest,
coordinates, cellvolumes, bounding_box, lower_corner, upper_corner,
voxel_centers, cells, scale`. Keep `add_shape`. Remove `World`.
Watch for name clashes with existing exports (`center`, `lengths`, `radius`,
`radii`, `half_height` are shape accessors; `center` is extended for `Region`,
the rest are not reused). `scale` shadows nothing in Base but check Aqua.

## 3. Conventions and constraints
- Julia 1.10 compat, Julia 1.12 locally. Run tests with
  `julia --project -e 'using Pkg; Pkg.test()'` from the repo root, or a single
  file with `julia --project -e 'using Test, VoxelShapes, StaticArrays; using LinearAlgebra: norm; include("test/test_xxx.jl")'`.
- CUDA tests are skipped when no GPU is present (check how `test/test_cuda.jl`
  guards); keep that behavior.
- Aqua tests run in `test/test_aqua.jl`; new exports must not introduce method
  ambiguities or piracy.
- Everything reachable from the GPU kernel must stay `isbits`: `Geometry` must be
  isbits when its shapes are (a `Tuple` of shapes, as `World` did). `Region` and
  `CompositeGrid` are host-only and are never passed to the kernel; only the
  converted `origin`/`voxel_size` tuples are.
- Docstrings for every public function, in the style of the existing code.
- Follow the existing module layout: one `module` per file under `src/`,
  `using ..Types`, exports re-exported from `src/VoxelShapes.jl`.
- Do not add dependencies. Do not add Gila as a dependency; the Gila-side
  extension (`GlaVol(::Region)`, `GlaCmpVol(::CompositeGrid)`) is out of scope
  for this repository.

## 4. Work breakdown

### Phase 1: Grids (Opus)
Files: new `src/Grids.jl`; `src/VoxelShapes.jl` (add `include("Grids.jl")`,
`using .Grids`, exports); new `test/test_grids.jl`; `test/runtests.jl`.
Deliver everything in section 2.1 **except** the shape overload of `refine`.
Tests: port every case in Gila's `test/cmpVolTest.jl` (fetch it from GitHub) to
the new names, including the 7-region reference case, the box-flush-against-a-face
case, nested refinement, anisotropic factors, and every invariant violation
(each must throw `ArgumentError`). Also test `coordinates` ordering against
`voxel_centers`, `cellvolumes` sums, `bounding_box`, `==`, `show`.

### Phase 2: Shape bounding boxes (Sonnet), parallel with Phase 1
Files: `src/Types.jl` (declare `bounding_box`, default all-infinite),
every shape module, `src/Transforms.jl`, `src/CSG.jl`, `src/VoxelShapes.jl`
(export `bounding_box` from the Types export line), `test/test_shapes.jl`,
`test/test_transforms.jl`, `test/test_csg.jl`.
Tests: for every shape, sample many random points; every point with
`point in shape` must lie inside the box (conservativeness), and for finite
shapes the box must be tight along each axis (touches the shape within a
tolerance) where the formula is exact (ellipsoid, cuboid, cylinder, torus,
capsule, cone, axis-aligned slab). Rotated cuboid by 45° about z must widen x/y
by `sqrt(2)`. CSG rules as specified.

### Phase 3: Geometry, rasterize, CompositeField, remove World (Opus)
Depends on Phases 1 and 2 being merged.
Files: new `src/Geometry.jl`, new `src/Rasterize.jl`, `src/VoxelShapes.jl`
(remove World, include the new files, exports), `src/Grids.jl` (add the
shape overload of `refine`), tests: rename `test/test_world.jl` to
`test/test_geometry.jl` and rewrite for `Geometry`/`Region`/`rasterize`; new
`test/test_composite.jl`; update `test/test_integration.jl`,
`test/test_cuda.jl`, `test/runtests.jl`.
Key tests:
- Single-region rasterize reproduces the old World results (translate the
  existing tests: 3×3×3 unit grid with sphere at (1.5,1.5,1.5) etc. Note the
  old tests used a corner origin of 0; the equivalent `Region` has center
  `(3//2,3//2,3//2)`).
- Composite: rasterize a sphere into the 7-region reference grid; `regrid` to
  the finest scale must equal `rasterize` of the same geometry into a single
  `Region` at that finest scale wherever the coarse regions are fully inside or
  fully outside the sphere, and the core block must equal it exactly.
- `vec(field)` length equals `length(grid)`; `regionview` shapes equal region
  cells; `eachregion` order matches `regions`.
- `refine(grid, shape; padding)` puts the shape's bounding box inside the core.
- `sample(geometry, pos, vs)` matches `rasterize` at voxel centers.
- CUDA path (guarded) for single region and composite.

### Phase 4: Docs, README, examples, version (Sonnet)
Depends on Phase 3.
Files: `Project.toml` (version `0.2.0`), `README.md`, `docs/make.jl` (add
`VoxelShapes.Grids` to `modules`, add a "Grids and refinement" page),
`docs/src/*.md` (every `World` reference: `index.md`, `gpu.md`, `api.md`,
`csg.md`, `antialiasing.md`, `transforms.md`; new `grids.md`; `extending.md`
gets `bounding_box`), `examples/01`–`07` (replace `World`/`Array(world)` with
`Geometry`/`Region`/`rasterize`; add `examples/08_refine.jl` showing a composite
grid around a small sphere, rendered via `regrid`). The examples need GLMakie,
which may not be installed; at minimum run the non-plotting part of each
example. Do not regenerate PNGs unless GLMakie is available. Build the docs
locally with `julia -e 'using Pkg; Pkg.activate("docs"); include("docs/make.jl")'`
if it works offline; otherwise check the markdown by eye.

### Phase 5: Verification (Sonnet)
Run the full test suite and Aqua; fix anything small; report anything large
back to the orchestrator instead of guessing.

## 5. Out of scope (later rounds)
- Gila-side package extension converting `Region`/`CompositeGrid` to
  `GlaVol`/`GlaCmpVol` and back.
- Making CUDA a package extension of VoxelShapes.
- Non-outward snapping modes.

## 6. Code clarity rules (apply while writing, audit before finishing)

Every phase must follow the "no hollow abstractions" rule. A named function is
a claim that the concept recurs or is complex enough that isolating it helps the
reader. When that claim is false the name is pure cost.

**Deletion test:** mentally inline the helper at its call sites. If the program
got easier to read, the helper was hollow. Remove it.

Three species to avoid:
1. **Shadow interface**: a private helper (`_foo`) that implements an operation
   whose canonical home is `Base.:(==)`, `Base.show`, `Base.vec`, `Base.length`,
   `Base.iterate`, etc., with the canonical method reduced to a one-line
   trampoline. Put the body in the canonical method; route other callers through
   the interface.
2. **One-shot fragment**: one caller, roughly five lines or fewer, underscore
   name. Inline it and keep a one-sentence comment above the block.
3. **Pass-through layer**: `f` only calls `g`, `g` has no other caller, one of
   them adds nothing. Collapse into one function under the better name.

Keep helpers that earn their name: dispatch points with methods on several
types, multiple genuine call sites, deep logic (roughly 15+ lines of gnarly
work), named domain concepts, and public or separately tested API.

Gila's `glaCmpVol.jl` has several tiny private helpers (`_lwrEdg`, `_uprEdg`,
`_ovrLap`, `_facTup`, `_boxBnd`, `_subVol`, `_snpBox`, `_crvVol`). Do **not**
port them one for one. `lower_corner`/`upper_corner` are public accessors with
many call sites and earn their names. `_ovrLap` is a two-liner used in two
places; judge it honestly. `_facTup`/`_boxBnd` are argument normalization and
belong inline in `refine` as `if`/`elseif` on the argument type, or as a small
dispatch set only if that reads better. `_snpBox`/`_crvVol`/`_subVol` are the
real algorithm; keep whichever of them are long enough to lower `refine`'s
cognitive load and inline the rest.

Also avoid: a validation helper that only wraps another validation helper, a
helper created just to hold one error check, and a helper extracted
preemptively "for cleanliness" before a second caller exists.

**Phase 5 (verification) must run this audit** over `src/Grids.jl`,
`src/Geometry.jl`, `src/Rasterize.jl`, and every shape file touched: list
private helpers, count call sites across the repo with grep, classify each,
inline or collapse the hollow ones, and re-run the tests (behavior must not
change). Report each decision with its call-site count.

## 7. Prose style (all phases) and the humanizer pass

All comments, docstrings, README and docs text must follow `STYLE.md` in the
repo root. It holds two samples of Paul's writing (thesis prose and Julia
docstrings), the habits to copy, and the AI tells to remove. Read it before
writing any prose. Not every function gets a docstring; follow the guidance
there on where detail belongs and where it is left out.

### Phase 6: Humanizer pass (Sonnet), after Phase 5
Read `/Users/pvirally/.claude/skills/humanizer/SKILL.md` and `STYLE.md`. Then
go through every file touched in this release (`src/*.jl`, `test/*.jl`
comments only, `README.md`, `docs/src/*.md`, `examples/*.jl` comments) and
rewrite comments, docstrings and documentation prose to match Paul's voice.
Rules: keep every fact and identifier, invent nothing, leave code untouched,
no em or en dashes anywhere in the final text, apply the docstring layout from
`STYLE.md` (signature line, one sentence, `# Arguments`/`# Fields`/`# Returns`
fragments without trailing periods). Delete comments that merely restate the
line below them. Re-run the tests after (doctests, if any, and `Pkg.test()`)
since docstrings can contain code. Report a short summary of the kinds of
changes made, not a diff.
