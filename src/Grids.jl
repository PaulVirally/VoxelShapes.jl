module Grids

export Region, CompositeGrid, refine
export cells, scale, lower_corner, upper_corner, voxel_centers
export regions, nregions, finest, coordinates, cellvolumes

using ..Types
import ..Types: bounding_box

"""
    Region

A solid, axis-aligned cuboid of uniform voxels.

A `Region` is the VoxelShapes counterpart of a `GlaVol` in
GilaElectromagnetics.jl: it stores a cell count, a voxel side length, and the
center of the cuboid, all in exact rational arithmetic so that regions of
different resolutions can be compared and glued together without round-off.
Lengths are unitless; when a grid is handed to Gila they are wavelengths.

The region spans `center - cells * scale / 2` to `center + cells * scale / 2`.
The center of cell `i` along an axis is `lower_corner + (i - 1/2) * scale`,
which is exactly the convention used by [`rasterize`](@ref VoxelShapes.rasterize) and by Gila.

# Fields
- `cells::NTuple{3,Int}`: number of voxels along each axis
- `scale::NTuple{3,Rational{Int}}`: voxel side lengths
- `center::NTuple{3,Rational{Int}}`: center of the region, `(0//1, 0//1, 0//1)` by default

Integer entries in `scale` and `center` are converted to `Rational{Int}`, so
`Region((2, 2, 2), (1, 1, 1))` is the same as `Region((2, 2, 2), (1//1, 1//1, 1//1))`.
An odd cell count is fine for a lone region. Only a [`CompositeGrid`](@ref)
demands even counts, the same way Gila does.
"""
struct Region
    cells::NTuple{3,Int}
    scale::NTuple{3,Rational{Int}}
    center::NTuple{3,Rational{Int}}

    function Region(cells::NTuple{3,Integer},
                    scale::NTuple{3,Union{Integer,Rational}},
                    center::NTuple{3,Union{Integer,Rational}}=(0 // 1, 0 // 1, 0 // 1))
        if any(cells .< 1)
            bad = findall(collect(cells) .< 1)
            throw(ArgumentError("A region needs at least one cell in every dimension, got $(cells) (dimension(s) $(bad) are empty)."))
        end
        if any(scale .<= 0)
            bad = findall(collect(scale) .<= 0)
            throw(ArgumentError("Voxel side lengths must be positive, got $(scale) (dimension(s) $(bad) are not)."))
        end
        return new(Int.(cells), Rational{Int}.(scale), Rational{Int}.(center))
    end
end

"""
    cells(region::Region) -> NTuple{3,Int}

The number of voxels of `region` along each axis.
"""
cells(region::Region) = region.cells

"""
    scale(region::Region) -> NTuple{3,Rational{Int}}

The voxel side lengths of `region` along each axis.
"""
scale(region::Region) = region.scale

Types.center(region::Region) = region.center

"""
    lower_corner(region::Region) -> NTuple{3,Rational{Int}}

The corner of `region` with the smallest coordinate along every axis, that is
`center(region) - cells(region) * scale(region) / 2`.
"""
lower_corner(region::Region) =
    ntuple(dir -> region.center[dir] - region.cells[dir] * region.scale[dir] // 2, 3)

"""
    upper_corner(region::Region) -> NTuple{3,Rational{Int}}

The corner of `region` with the largest coordinate along every axis, that is
`center(region) + cells(region) * scale(region) / 2`.
"""
upper_corner(region::Region) =
    ntuple(dir -> region.center[dir] + region.cells[dir] * region.scale[dir] // 2, 3)

"""
    voxel_centers(region::Region) -> NTuple{3,StepRange{Rational{Int},Rational{Int}}}

The voxel center coordinates of `region`, one exact range per axis.

Range `dir` has `cells(region)[dir]` entries, starts at
`lower_corner(region)[dir] + scale(region)[dir] / 2`, and steps by
`scale(region)[dir]`. These are the same values as the `grd` field of a Gila
`GlaVol`.
"""
function voxel_centers(region::Region)
    lwr = lower_corner(region)
    return ntuple(3) do dir
        stp = region.scale[dir]
        fst = lwr[dir] + stp // 2
        StepRange(fst, stp, fst + (region.cells[dir] - 1) * stp)
    end
end

Base.size(region::Region) = region.cells
Base.length(region::Region) = prod(region.cells)

Base.:(==)(regA::Region, regB::Region) =
    regA.cells == regB.cells && regA.scale == regB.scale && regA.center == regB.center

Base.hash(region::Region, h::UInt) =
    hash(region.center, hash(region.scale, hash(region.cells, hash(:Region, h))))

function Base.show(io::IO, region::Region)
    print(io, "Region(", join(region.cells, "×"), " cells, scale (",
          join(region.scale, ", "), "), center (", join(region.center, ", "), "))")
end
Base.show(io::IO, ::MIME"text/plain", region::Region) = show(io, region)

"""
    CompositeGrid

A set of disjoint [`Region`](@ref)s, possibly of different resolutions, that
exactly tile one rectangular domain but behave like a single grid.

The order of the regions is meaningful: it fixes the flat layout used for
fields over the grid (region blocks concatenated in `regions` order, column
major inside a block) and matches the degree of freedom layout of Gila's
`GlaCmpVol`.

# Fields
- `regions::Vector{Region}`: the regions of the tiling, in layout order

# Invariants

The inner constructor checks, and throws `ArgumentError` if any of them fails,
that the regions

1. are nonempty in number and have at least one cell per dimension,
2. have even cell counts in every dimension,
3. have pairwise commensurate scales (one an integer multiple of the other),
4. have disjoint interiors (face, edge, and corner contact is allowed),
5. sit on a common grid (corners differ by a whole number of common cells),
6. satisfy the pairwise partition parity condition Gila's external
   interaction bookkeeping needs, and
7. exactly tile the bounding box of their union, with no gaps.
"""
struct CompositeGrid
    regions::Vector{Region}

    function CompositeGrid(regs::Vector{Region})
        _check_grid(regs)
        return new(copy(regs))
    end
end

"""
    CompositeGrid(region::Region)

Construct the one-region composite grid around a single [`Region`](@ref).
"""
CompositeGrid(region::Region) = CompositeGrid([region])

"""
    regions(grid::CompositeGrid) -> Vector{Region}

The regions of `grid`, in flat layout order.
"""
regions(grid::CompositeGrid) = grid.regions

"""
    nregions(grid::CompositeGrid) -> Int

The number of regions in `grid`.
"""
nregions(grid::CompositeGrid) = length(grid.regions)

"""
    finest(grid::CompositeGrid) -> NTuple{3,Rational{Int}}

The smallest voxel side length present in `grid`, taken elementwise over the
scales of its regions.
"""
finest(grid::CompositeGrid) =
    reduce((sclA, sclB) -> min.(sclA, sclB), (reg.scale for reg in grid.regions))

"""
    bounding_box(grid::CompositeGrid) -> (lower, upper)

The corners of the cuboid tiled by `grid`, as a pair of
`NTuple{3,Rational{Int}}`.
"""
function bounding_box(grid::CompositeGrid)
    lwr = reduce((a, b) -> min.(a, b), lower_corner.(grid.regions))
    upr = reduce((a, b) -> max.(a, b), upper_corner.(grid.regions))
    return (Tuple(lwr), Tuple(upr))
end

Base.length(grid::CompositeGrid) = sum(length, grid.regions)

Base.:(==)(gridA::CompositeGrid, gridB::CompositeGrid) = gridA.regions == gridB.regions

function Base.show(io::IO, grid::CompositeGrid)
    num = nregions(grid)
    print(io, "CompositeGrid ($num region", num == 1 ? "" : "s", ")")
    for (idx, reg) in enumerate(grid.regions)
        print(io, "\n  [$idx] (", join(reg.cells, "×"), ") cells, scale (",
              join(reg.scale, ", "), "), center (", join(reg.center, ", "), ")")
    end
end
Base.show(io::IO, ::MIME"text/plain", grid::CompositeGrid) = show(io, grid)

function _check_grid(regs::Vector{Region})
    if isempty(regs)
        throw(ArgumentError("A composite grid needs at least one region."))
    end
    for (idx, reg) in enumerate(regs)
        if any(isodd.(reg.cells))
            bad = findall(collect(isodd.(reg.cells)))
            throw(ArgumentError("Region $idx has an odd number of cells in dimension(s) $(bad): $(reg.cells). Every region of a composite grid must have an even cell count in every dimension."))
        end
    end
    gcdScl = reduce((sclA, sclB) -> gcd.(sclA, sclB), (reg.scale for reg in regs))
    for idxA in 1:length(regs), idxB in (idxA + 1):length(regs)
        regA, regB = regs[idxA], regs[idxB]
        _check_common_scale(regA, regB, idxA, idxB)
        # Interiors share volume. Contact along a face, an edge, or a corner is
        # fine: touching regions are exactly what a refined tiling is made of.
        lwr = max.(lower_corner(regA), lower_corner(regB))
        upr = min.(upper_corner(regA), upper_corner(regB))
        if all(lwr .< upr)
            throw(ArgumentError("Regions $idxA and $idxB overlap: region $idxA spans $(lower_corner(regA)) to $(upper_corner(regA)) and region $idxB spans $(lower_corner(regB)) to $(upper_corner(regB)). Regions may touch, but their interiors must be disjoint."))
        end
        offset = (lower_corner(regA) .- lower_corner(regB)) .// gcdScl
        if any(.!isinteger.(offset))
            throw(ArgumentError("Regions $idxA and $idxB do not share a common grid: their lower corners differ by $(lower_corner(regA) .- lower_corner(regB)), which is $(offset) common cells of size $(gcdScl). Region corners must lie on the common grid."))
        end
        _check_parity(regA, regB, idxA, idxB)
    end
    _check_tiling(regs)
    return nothing
end

#= Larger over smaller has to be a whole number in every dimension, otherwise no
common grid exists for the pair. =#
function _check_common_scale(regA::Region, regB::Region, idxA::Integer, idxB::Integer)
    ratio = max.(regA.scale, regB.scale) .// min.(regA.scale, regB.scale)
    if any(.!isinteger.(ratio))
        bad = findall(collect(.!isinteger.(ratio)))
        throw(ArgumentError("Regions $idxA and $idxB have incommensurate voxel scales in dimension(s) $(bad): $(regA.scale) and $(regB.scale) give the ratio(s) $(ratio). One scale must be an integer multiple of the other in every dimension."))
    end
    return nothing
end

#= Count the cells of each region in units of the coarsest common cell and
require the sum to be even, the condition Gila's external interaction
bookkeeping imposes on every pair of regions. =#
function _check_parity(regA::Region, regB::Region, idxA::Integer, idxB::Integer)
    maxScl = lcm.(regA.scale, regB.scale)
    divA = ntuple(dir -> maxScl[dir] ÷ regA.scale[dir], 3)
    divB = ntuple(dir -> maxScl[dir] ÷ regB.scale[dir], 3)
    for (idx, reg, divPar) in ((idxA, regA, divA), (idxB, regB, divB))
        if any(mod.(reg.cells, divPar) .!= 0)
            bad = findall(collect(mod.(reg.cells, divPar) .!= 0))
            other = idx == idxA ? idxB : idxA
            throw(ArgumentError("Region $idx cannot be partitioned against region $other in dimension(s) $(bad): $(reg.cells) cells do not divide into $(divPar) sub-lattices. Adjust the refinement box so that the region spans a whole number of coarse cells."))
        end
    end
    parA = regA.cells .÷ divA
    parB = regB.cells .÷ divB
    total = parA .+ parB
    if any(isodd.(total))
        bad = findall(collect(isodd.(total)))
        throw(ArgumentError("Regions $idxA and $idxB violate the partition parity condition in dimension(s) $(bad): $(parA) cells per partition in region $idxA plus $(parB) in region $idxB gives $(total). Move the refinement box or change the refinement factor so that every dimension sums to an even number."))
    end
    return nothing
end

#= Disjoint regions whose volumes add up to the volume of their bounding box
tile that box exactly. =#
function _check_tiling(regs::Vector{Region})
    minEdg = reduce((a, b) -> min.(a, b), lower_corner.(regs))
    maxEdg = reduce((a, b) -> max.(a, b), upper_corner.(regs))
    bndVol = prod(maxEdg .- minEdg)
    regVol = sum(prod(reg.cells .* reg.scale) for reg in regs)
    if regVol != bndVol
        throw(ArgumentError("The regions do not tile their bounding box: they fill $(regVol) of a box of $(bndVol) spanning $(Tuple(minEdg)) to $(Tuple(maxEdg)). A composite grid is one solid cuboid, so gaps are not allowed. Use one composite grid per body."))
    end
    return nothing
end

_factor_tuple(factor::Integer) = ntuple(_ -> Int(factor), 3)
_factor_tuple(factor::NTuple{3,Integer}) = Int.(factor)
_factor_tuple(factor) =
    throw(ArgumentError("A refinement factor must be an Integer or an NTuple{3,Integer}, got $(typeof(factor))."))

# Corners of the box to refine inside
_box_bounds(box::Region) = (lower_corner(box), upper_corner(box))
_box_bounds(box::Tuple{NTuple{3,Real},NTuple{3,Real}}) =
    (box[1] .- box[2] ./ 2, box[1] .+ box[2] ./ 2)
_box_bounds(box) =
    throw(ArgumentError("A refinement box must be a Region or a tuple (center, dims) of center and full side lengths, got $(typeof(box))."))

#= Cell index bounds of a box inside a region, half open and zero based. The
bounds grow outward to cell boundaries and then to even cell indices, keeping
the core and every complement slab at an even cell count. =#
function _snap_box(region::Region, boxLwr, boxUpr)
    regLwr = lower_corner(region)
    rawLwr = ntuple(dir -> floor(Int, (boxLwr[dir] - regLwr[dir]) / region.scale[dir]), 3)
    rawUpr = ntuple(dir -> ceil(Int, (boxUpr[dir] - regLwr[dir]) / region.scale[dir]), 3)
    idxLwr = ntuple(dir -> 2 * fld(max(rawLwr[dir], 0), 2), 3)
    idxUpr = ntuple(dir -> 2 * cld(min(rawUpr[dir], region.cells[dir]), 2), 3)
    return idxLwr, idxUpr
end

#= The piece of region spanning cells idxLwr+1 through idxUpr, at the parent
voxel scale divided by fac. =#
function _sub_region(region::Region, idxLwr::NTuple{3,Int}, idxUpr::NTuple{3,Int},
                     fac::NTuple{3,Int}=(1, 1, 1))
    regLwr = lower_corner(region)
    lwr = ntuple(dir -> regLwr[dir] + idxLwr[dir] * region.scale[dir], 3)
    upr = ntuple(dir -> regLwr[dir] + idxUpr[dir] * region.scale[dir], 3)
    cel = ntuple(dir -> (idxUpr[dir] - idxLwr[dir]) * fac[dir], 3)
    scl = ntuple(dir -> region.scale[dir] // fac[dir], 3)
    ctr = ntuple(dir -> (lwr[dir] + upr[dir]) // 2, 3)
    return Region(cel, scl, ctr)
end

#= Split a region into the refined core plus the complement slabs, in the fixed
order core, xlo, xhi, ylo, yhi, zlo, zhi. The x slabs span the full cross
section, the y slabs only the middle in x, and the z slabs only the middle in x
and y, so the pieces are disjoint and fill the region. Empty slabs are
dropped. =#
function _carve(region::Region, boxLwr, boxUpr, fac::NTuple{3,Int})
    idxLwr, idxUpr = _snap_box(region, boxLwr, boxUpr)
    cel = region.cells
    pcs = Region[_sub_region(region, idxLwr, idxUpr, fac)]
    if idxLwr[1] > 0
        push!(pcs, _sub_region(region, (0, 0, 0), (idxLwr[1], cel[2], cel[3])))
    end
    if idxUpr[1] < cel[1]
        push!(pcs, _sub_region(region, (idxUpr[1], 0, 0), (cel[1], cel[2], cel[3])))
    end
    if idxLwr[2] > 0
        push!(pcs, _sub_region(region, (idxLwr[1], 0, 0), (idxUpr[1], idxLwr[2], cel[3])))
    end
    if idxUpr[2] < cel[2]
        push!(pcs, _sub_region(region, (idxLwr[1], idxUpr[2], 0), (idxUpr[1], cel[2], cel[3])))
    end
    if idxLwr[3] > 0
        push!(pcs, _sub_region(region, (idxLwr[1], idxLwr[2], 0), (idxUpr[1], idxUpr[2], idxLwr[3])))
    end
    if idxUpr[3] < cel[3]
        push!(pcs, _sub_region(region, (idxLwr[1], idxLwr[2], idxUpr[3]), (idxUpr[1], idxUpr[2], cel[3])))
    end
    return pcs
end

"""
    refine(grid, box; factor=2, snap=:outward) -> CompositeGrid

Refine a grid inside a box by carving.

Every region whose interior meets `box` is replaced by a refined core plus the
complement slabs that fill the rest of the region. The core takes `factor` times
as many voxels at `factor` times the resolution over the same physical extent.
Regions the box misses are left alone and keep their position in the region
list, so the pieces of region `i` appear where region `i` used to be: the core
first, then the nonempty slabs in the order xlo, xhi, ylo, yhi, zlo, zhi. This
order fixes the flat field layout, and it is identical to the one Gila's
`refine` produces.

Snapping happens per region. The box is first intersected with the region, then
grown outward to that region's voxel boundaries, then grown outward again to
even cell indices so that the core and every slab keep an even cell count. A box
that crosses regions of different scale therefore snaps differently in each of
them. The core always covers at least the requested box.

Refinement may be applied repeatedly, including in a nested fashion.

# Arguments
- `grid`: the [`CompositeGrid`](@ref) or [`Region`](@ref) to refine; a `Region`
  is wrapped in a one-region composite grid first
- `box`: the refinement box, either a [`Region`](@ref) whose outer bounds are
  used, or a tuple `(center, dims)` of center and full side lengths
- `factor=2`: the refinement factor, an `Integer` or an `NTuple{3,Integer}` with
  every entry at least 1. An entry of 1 leaves that direction at the parent
  resolution while still carving the region.
- `snap=:outward`: the snapping mode. Only `:outward` exists.

# Returns
- `CompositeGrid`: the refined grid, validated by the `CompositeGrid`
  constructor

# Throws
- `ArgumentError`: if `snap` is not `:outward`, if any refinement factor is less
  than 1, if `box` is not a supported type, or if the resulting tiling violates
  an invariant of [`CompositeGrid`](@ref)
"""
function refine end

function refine(grid::CompositeGrid, box; factor=2, snap::Symbol=:outward)
    if snap !== :outward
        throw(ArgumentError("Unknown snapping mode :$(snap). Only :outward is implemented."))
    end
    fac = _factor_tuple(factor)
    if any(fac .< 1)
        throw(ArgumentError("Refinement factors must be at least 1, got $(fac)."))
    end
    boxLwr, boxUpr = _box_bounds(box)
    return _refine_corners(grid, boxLwr, boxUpr, fac)
end

refine(region::Region, box; factor=2, snap::Symbol=:outward) =
    refine(CompositeGrid(region), box; factor=factor, snap=snap)

"""
    refine(grid, shape::AbstractFillableShape; factor=2, padding=0, snap=:outward)

Refine a grid around a shape.

The refinement box is `bounding_box(shape)` grown by `padding` on every side and
then clamped to the bounding box of `grid`. A shape that does not reach into the
grid at all, or one that lies entirely outside it, leaves the grid unchanged.
Shapes that do not implement `bounding_box` report an infinite box, so the whole
grid is refined.

Note that the padding is a physical length, not a number of voxels. A padding of
one coarse voxel is a reasonable choice when the shape has an anti-aliased
boundary, since the boundary voxels are the ones that gain from the extra
resolution.

# Arguments
- `grid`: the [`CompositeGrid`](@ref) or [`Region`](@ref) to refine
- `shape`: the shape whose bounding box is refined
- `factor=2`: the refinement factor, an `Integer` or an `NTuple{3,Integer}`
- `padding=0`: extra length added on every side of the bounding box, a `Real` or
  an `NTuple{3,Real}` giving one padding per axis
- `snap=:outward`: the snapping mode. Only `:outward` exists.
"""
function refine(grid::CompositeGrid, shape::AbstractFillableShape;
                factor=2, padding=0, snap::Symbol=:outward)
    if snap !== :outward
        throw(ArgumentError("Unknown snapping mode :$(snap). Only :outward is implemented."))
    end
    fac = _factor_tuple(factor)
    if any(fac .< 1)
        throw(ArgumentError("Refinement factors must be at least 1, got $(fac)."))
    end
    pad = padding isa Real ? ntuple(_ -> padding, 3) : padding
    shpLwr, shpUpr = bounding_box(shape)
    grdLwr, grdUpr = bounding_box(grid)
    boxLwr = max.(shpLwr .- pad, grdLwr)
    boxUpr = min.(shpUpr .+ pad, grdUpr)
    all(boxLwr .< boxUpr) || return grid
    return _refine_corners(grid, boxLwr, boxUpr, fac)
end

refine(region::Region, shape::AbstractFillableShape; factor=2, padding=0, snap::Symbol=:outward) =
    refine(CompositeGrid(region), shape; factor=factor, padding=padding, snap=snap)

#= Carve every region the box reaches into, and leave the rest in place. The
box corners may be Rational or Float64; _snap_box only needs floor and ceil. =#
function _refine_corners(grid::CompositeGrid, boxLwr, boxUpr, fac::NTuple{3,Int})
    newRegs = Region[]
    for reg in grid.regions
        lwr = max.(boxLwr, lower_corner(reg))
        upr = min.(boxUpr, upper_corner(reg))
        if all(lwr .< upr)
            append!(newRegs, _carve(reg, lwr, upr, fac))
        else
            push!(newRegs, reg)
        end
    end
    return CompositeGrid(newRegs)
end

struct CompositeCoordinates
    grid::CompositeGrid
end

"""
    coordinates(grid::CompositeGrid)

Iterate over the voxels of `grid` in flat layout order.

The iterator yields one `(pos, cellvolume, regionindex)` triple per voxel:
`pos::NTuple{3,Float64}` is the voxel center, `cellvolume::Float64` is the
volume of that voxel, and `regionindex::Int` is the index of the region the
voxel belongs to. Voxels come out region by region and, inside a region, in
column major order over the voxel grid, which is the order a field over the grid
is stored in.
"""
coordinates(grid::CompositeGrid) = CompositeCoordinates(grid)

Base.eltype(::Type{CompositeCoordinates}) = Tuple{NTuple{3,Float64},Float64,Int}
Base.length(itr::CompositeCoordinates) = length(itr.grid)
Base.IteratorSize(::Type{CompositeCoordinates}) = Base.HasLength()

function Base.iterate(itr::CompositeCoordinates, stt::Tuple{Int,Int}=(1, 1))
    regIdx, celIdx = stt
    regs = itr.grid.regions
    regIdx > length(regs) && return nothing
    reg = regs[regIdx]
    if celIdx > length(reg)
        return iterate(itr, (regIdx + 1, 1))
    end
    lwr = lower_corner(reg)
    crtInd = CartesianIndices(reg.cells)[celIdx]
    pos = ntuple(dir -> Float64(lwr[dir] + (crtInd[dir] - 1 // 2) * reg.scale[dir]), 3)
    return ((pos, Float64(prod(reg.scale)), regIdx), (regIdx, celIdx + 1))
end

"""
    cellvolumes(grid::CompositeGrid) -> Vector{Float64}

The voxel volumes of `grid`, one entry per voxel, in flat layout order.

The result has length `length(grid)` and lines up entry by entry with a scalar
field over the grid. A region has a single voxel size, so the whole block of a
region carries one value.

Note that this differs from Gila's `cellvolumes`, which repeats each volume
three times because Gila's fields are vector valued. VoxelShapes fields are
scalar, so there is exactly one entry per voxel here.
"""
function cellvolumes(grid::CompositeGrid)
    vols = Float64[]
    for reg in grid.regions
        append!(vols, fill(Float64(prod(reg.scale)), length(reg)))
    end
    return vols
end

end # module Grids
