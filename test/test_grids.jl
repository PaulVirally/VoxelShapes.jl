# Grids: Region and CompositeGrid geometry.

const grdScl16 = (1 // 16, 1 // 16, 1 // 16)
const grdScl32 = (1 // 32, 1 // 32, 1 // 32)
const grdCtr0 = (0 // 1, 0 // 1, 0 // 1)

# Base region for the carving tests: 8 cells of 1/16 per side, centered on the
# origin, so it spans -1/4 to 1/4 in every dimension.
mkgrid() = Region((8, 8, 8), grdScl16, grdCtr0)

# Exact rational volume filled by the regions of a composite grid
tiled_volume(g) = sum(prod(cells(r) .* scale(r)) for r in regions(g))

# True if no two regions of a composite grid share interior
function disjoint_regions(g)
    regs = regions(g)
    for i in 1:length(regs), j in (i + 1):length(regs)
        lwr = max.(lower_corner(regs[i]), lower_corner(regs[j]))
        upr = min.(upper_corner(regs[i]), upper_corner(regs[j]))
        all(lwr .< upr) && return false
    end
    return true
end

@testset "Region" begin
    r = mkgrid()
    @test cells(r) == (8, 8, 8)
    @test scale(r) == grdScl16
    @test center(r) == grdCtr0
    @test size(r) == (8, 8, 8)
    @test length(r) == 512
    @test lower_corner(r) == (-1 // 4, -1 // 4, -1 // 4)
    @test upper_corner(r) == (1 // 4, 1 // 4, 1 // 4)

    # Integer entries in scale and center are promoted to Rational{Int}
    ri = Region((2, 2, 2), (1, 1, 1), (0, 0, 0))
    @test scale(ri) === (1 // 1, 1 // 1, 1 // 1)
    @test center(ri) === (0 // 1, 0 // 1, 0 // 1)
    @test lower_corner(ri) == (-1 // 1, -1 // 1, -1 // 1)

    # Default center is the origin
    @test Region((2, 2, 2), grdScl16) == Region((2, 2, 2), grdScl16, grdCtr0)

    # Voxel centers: cells(r) values per axis, offset half a voxel from the
    # corner
    vc = voxel_centers(r)
    @test length(vc) == 3
    for dir in 1:3
        @test length(vc[dir]) == cells(r)[dir]
        @test first(vc[dir]) == lower_corner(r)[dir] + scale(r)[dir] // 2
        @test last(vc[dir]) == upper_corner(r)[dir] - scale(r)[dir] // 2
        @test step(vc[dir]) == scale(r)[dir]
        @test eltype(vc[dir]) == Rational{Int}
    end
    # An odd cell count is fine for a lone region
    rodd = Region((3, 4, 5), grdScl16, grdCtr0)
    @test cells(rodd) == (3, 4, 5)
    @test length(voxel_centers(rodd)[1]) == 3
    @test voxel_centers(rodd)[1][2] == 0 // 1

    # Anisotropic scales
    ra = Region((2, 4, 6), (1 // 2, 1 // 4, 1 // 8), (1 // 1, 0 // 1, -1 // 1))
    @test lower_corner(ra) == (1 // 1 - 1 // 2, -1 // 2, -1 // 1 - 3 // 8)
    @test upper_corner(ra) == (1 // 1 + 1 // 2, 1 // 2, -1 // 1 + 3 // 8)

    # Equality and show
    @test mkgrid() == mkgrid()
    @test mkgrid() != Region((8, 8, 8), grdScl32, grdCtr0)
    @test mkgrid() != Region((8, 8, 8), grdScl16, (1 // 2, 0 // 1, 0 // 1))
    @test hash(mkgrid()) == hash(mkgrid())
    str = sprint(show, r)
    @test occursin("Region", str)
    @test occursin("8×8×8", str)
    @test occursin("1//16", str)
    @test sprint(show, MIME"text/plain"(), r) == str

    # Rejections
    @test_throws ArgumentError Region((0, 4, 4), grdScl16, grdCtr0)
    @test_throws ArgumentError Region((4, -1, 4), grdScl16, grdCtr0)
    @test_throws ArgumentError Region((4, 4, 4), (0 // 1, 1 // 16, 1 // 16), grdCtr0)
    @test_throws ArgumentError Region((4, 4, 4), (-1 // 16, 1 // 16, 1 // 16), grdCtr0)
end

@testset "CompositeGrid trivial composite" begin
    r = mkgrid()
    g = CompositeGrid(r)
    @test nregions(g) == 1
    @test regions(g) == [r]
    @test g == CompositeGrid([r])
    @test finest(g) == grdScl16
    @test length(g) == prod(cells(r))
    @test bounding_box(g) == (lower_corner(r), upper_corner(r))

    # One entry per cell
    crd = collect(coordinates(g))
    @test length(crd) == prod(cells(r))
    @test length(coordinates(g)) == prod(cells(r))
    @test eltype(coordinates(g)) == Tuple{NTuple{3,Float64},Float64,Int}
    @test Base.IteratorSize(typeof(coordinates(g))) == Base.HasLength()
    @test length(cellvolumes(g)) == prod(cells(r))
    @test all(trp -> trp[3] == 1, crd)
    domVol = Float64(prod(cells(r) .* scale(r)))
    @test sum(trp -> trp[2], crd) ≈ domVol
    @test sum(cellvolumes(g)) ≈ domVol

    # Every voxel center of the region shows up exactly once
    vc = voxel_centers(r)
    @test Set(first.(crd)) ==
          Set(ntuple(dir -> Float64(vc[dir][ind[dir]]), 3)
              for ind in CartesianIndices(cells(r)))

    # Show
    str = sprint(show, g)
    @test occursin("CompositeGrid (1 region)", str)
    @test occursin("[1]", str)
    @test sprint(show, MIME"text/plain"(), g) == str
end

@testset "CompositeGrid single refine" begin
    r = mkgrid()
    g = refine(CompositeGrid(r), (grdCtr0, (1 // 8, 1 // 8, 1 // 8)))
    # Core plus six slabs
    @test nregions(g) == 7
    @test disjoint_regions(g)
    @test tiled_volume(g) == prod(cells(r) .* scale(r))
    core = regions(g)[1]
    @test scale(core) == grdScl32
    @test cells(core) == (8, 8, 8)
    @test all(mod.(cells(core), 4) .== 0)
    @test center(core) == grdCtr0
    @test finest(g) == grdScl32
    # The slabs stay at the parent scale and keep even cell counts
    for slb in regions(g)[2:end]
        @test scale(slb) == grdScl16
        @test all(iseven.(cells(slb)))
        @test all(cells(slb) .> 0)
    end
    # Fixed slab order: xlo, xhi, ylo, yhi, zlo, zhi
    @test [cells(reg) for reg in regions(g)[2:end]] ==
          [(2, 8, 8), (2, 8, 8), (4, 2, 8), (4, 2, 8), (4, 4, 2), (4, 4, 2)]
    # The bounding box of the tiling is the original region
    @test bounding_box(g) == (lower_corner(r), upper_corner(r))
    # Cell count is conserved by the refinement bookkeeping
    @test length(cellvolumes(g)) == sum(prod(cells(reg)) for reg in regions(g))
    @test length(g) == sum(prod(cells(reg)) for reg in regions(g))
    @test occursin("CompositeGrid (7 regions)", sprint(show, g))
end

@testset "CompositeGrid box flush against a face" begin
    r = mkgrid()
    # Box hugs the low x face of the domain, so the x-low slab is empty
    g = refine(CompositeGrid(r), ((-3 // 16, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8)))
    @test nregions(g) == 6
    @test disjoint_regions(g)
    @test tiled_volume(g) == prod(cells(r) .* scale(r))
    core = regions(g)[1]
    @test scale(core) == grdScl32
    @test lower_corner(core)[1] == lower_corner(r)[1]
    @test all(reg -> all(cells(reg) .> 0), regions(g))
end

@testset "CompositeGrid box covering the domain" begin
    r = mkgrid()
    g = refine(CompositeGrid(r), (grdCtr0, (1 // 1, 1 // 1, 1 // 1)))
    @test nregions(g) == 1
    core = regions(g)[1]
    @test cells(core) == (16, 16, 16)
    @test scale(core) == grdScl32
    @test center(core) == grdCtr0
    @test tiled_volume(g) == prod(cells(r) .* scale(r))
    # A box that misses everything leaves the tiling alone
    away = refine(CompositeGrid(r), ((10 // 1, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8)))
    @test away == CompositeGrid(r)
end

@testset "CompositeGrid outward snapping" begin
    r = mkgrid()
    # Box spans 0 to 1/32 per dimension: off the cell grid and an odd cell count
    boxCtr = (1 // 64, 1 // 64, 1 // 64)
    boxDim = (1 // 32, 1 // 32, 1 // 32)
    g = refine(CompositeGrid(r), (boxCtr, boxDim))
    core = regions(g)[1]
    boxLwr = boxCtr .- boxDim .// 2
    boxUpr = boxCtr .+ boxDim .// 2
    # The core covers at least the requested box
    @test all(lower_corner(core) .<= boxLwr)
    @test all(upper_corner(core) .>= boxUpr)
    # Two parent cells per dimension, which is even, so every slab is even too
    @test all(cells(core) .÷ 2 .== 2)
    @test all(iseven.(cells(core) .÷ 2))
    @test all(reg -> all(iseven.(cells(reg))), regions(g))
    @test tiled_volume(g) == prod(cells(r) .* scale(r))
    @test disjoint_regions(g)
    # Snapping is per region, so the core sits on the parent grid
    @test all(isinteger.((lower_corner(core) .- lower_corner(r)) .// scale(r)))
end

@testset "CompositeGrid chained refine" begin
    r = mkgrid()
    g = refine(CompositeGrid(r), (grdCtr0, (1 // 8, 1 // 8, 1 // 8)))
    # The core spans -1/8 to 1/8; refine its high octant again
    g2 = refine(g, ((1 // 16, 1 // 16, 1 // 16), (1 // 8, 1 // 8, 1 // 8)))
    @test nregions(g2) == 10
    @test disjoint_regions(g2)
    @test tiled_volume(g2) == prod(cells(r) .* scale(r))
    @test finest(g2) == (1 // 64, 1 // 64, 1 // 64)
    inner = regions(g2)[1]
    @test scale(inner) == (1 // 64, 1 // 64, 1 // 64)
    @test cells(inner) == (8, 8, 8)
    # Three levels of resolution coexist, and the untouched slabs keep their
    # position at the end of the list
    @test Set(scale(reg) for reg in regions(g2)) ==
          Set([(1 // 64, 1 // 64, 1 // 64), grdScl32, grdScl16])
    @test [cells(reg) for reg in regions(g2)[5:end]] ==
          [cells(reg) for reg in regions(g)[2:end]]
    @test bounding_box(g2) == bounding_box(g)

    # A sub-box centered on the core leaves 1/32 slabs two cells wide. That
    # gives an odd cell count per partition against the 1/16 slabs.
    @test_throws ArgumentError refine(g, (grdCtr0, (1 // 16, 1 // 16, 1 // 16)))
end

@testset "CompositeGrid anisotropic factor" begin
    r = mkgrid()
    box = (grdCtr0, (1 // 8, 1 // 8, 1 // 8))
    gz = refine(CompositeGrid(r), box; factor=(1, 1, 2))
    coreZ = regions(gz)[1]
    @test scale(coreZ) == (1 // 16, 1 // 16, 1 // 32)
    @test cells(coreZ) == (4, 4, 8)
    @test nregions(gz) == 7
    @test tiled_volume(gz) == prod(cells(r) .* scale(r))
    @test disjoint_regions(gz)
    @test finest(gz) == (1 // 16, 1 // 16, 1 // 32)

    ga = refine(CompositeGrid(r), box; factor=(2, 2, 4))
    coreA = regions(ga)[1]
    @test scale(coreA) == (1 // 32, 1 // 32, 1 // 64)
    @test cells(coreA) == (8, 8, 16)
    @test tiled_volume(ga) == prod(cells(r) .* scale(r))
    @test disjoint_regions(ga)

    # A factor of 1 everywhere still carves the region
    gone = refine(CompositeGrid(r), box; factor=1)
    @test nregions(gone) == 7
    @test all(reg -> scale(reg) == grdScl16, regions(gone))
    @test tiled_volume(gone) == prod(cells(r) .* scale(r))
end

@testset "CompositeGrid refine argument errors" begin
    r = mkgrid()
    box = (grdCtr0, (1 // 8, 1 // 8, 1 // 8))
    @test_throws ArgumentError refine(CompositeGrid(r), box; snap=:inward)
    @test_throws ArgumentError refine(CompositeGrid(r), box; factor=0)
    @test_throws ArgumentError refine(CompositeGrid(r), box; factor=(2, 0, 2))
    @test_throws ArgumentError refine(CompositeGrid(r), box; factor=2.0)
    @test_throws ArgumentError refine(CompositeGrid(r), (0.0, 1.0))
    @test_throws ArgumentError refine(CompositeGrid(r), "everywhere")
end

@testset "CompositeGrid constructor rejections" begin
    # No regions at all
    @test_throws ArgumentError CompositeGrid(Region[])
    # Odd cell count
    @test_throws ArgumentError CompositeGrid(Region((3, 4, 4), grdScl32, grdCtr0))
    # Overlapping interiors
    rA = Region((4, 4, 4), grdScl32, grdCtr0)
    @test_throws ArgumentError CompositeGrid([rA, rA])
    # Incommensurate scales: 1/16 against 1/24
    rBad = Region((4, 4, 4), (1 // 24, 1 // 24, 1 // 24), (1 // 1, 0 // 1, 0 // 1))
    @test_throws ArgumentError CompositeGrid([Region((4, 4, 4), grdScl16, grdCtr0), rBad])
    # Misaligned grids: neighbour shifted by half a cell
    rOff = Region((4, 4, 4), grdScl32, (9 // 64, 0 // 1, 0 // 1))
    @test_throws ArgumentError CompositeGrid([rA, rOff])
    # Aligned and disjoint, but with a gap, so the regions do not tile
    rGap = Region((4, 4, 4), grdScl32, (1 // 4, 0 // 1, 0 // 1))
    @test_throws ArgumentError CompositeGrid([rA, rGap])
    # Partition parity violation: a (2,2,2) coarse region against a (2,2,2)
    # region at half the scale gives one cell per partition on both sides
    rCrs = Region((2, 2, 2), grdScl16, grdCtr0)
    rFin = Region((2, 2, 2), grdScl32, (1 // 16 + 1 // 32, 0 // 1, 0 // 1))
    @test_throws ArgumentError CompositeGrid([rCrs, rFin])
end

@testset "CompositeGrid flat layout" begin
    # Two touching regions of different scale, both 1/8 per side
    rCrs = Region((4, 4, 4), grdScl32, grdCtr0)
    rFin = Region((8, 8, 8), (1 // 64, 1 // 64, 1 // 64), (1 // 8, 0 // 1, 0 // 1))
    g = CompositeGrid([rCrs, rFin])
    @test nregions(g) == 2
    @test tiled_volume(g) == 2 * prod(cells(rCrs) .* scale(rCrs))
    @test finest(g) == (1 // 64, 1 // 64, 1 // 64)
    @test bounding_box(g) == ((-1 // 16, -1 // 16, -1 // 16), (3 // 16, 1 // 16, 1 // 16))
    @test length(g) == prod(cells(rCrs)) + prod(cells(rFin))

    crd = collect(coordinates(g))
    @test length(crd) == prod(cells(rCrs)) + prod(cells(rFin))
    # The first block is region 1, in column major order over its own grid
    fstBlk = crd[1:prod(cells(rCrs))]
    sndBlk = crd[(prod(cells(rCrs)) + 1):end]
    @test all(trp -> trp[3] == 1, fstBlk)
    @test all(trp -> trp[3] == 2, sndBlk)
    vcCrs = voxel_centers(rCrs)
    expPos = [ntuple(dir -> Float64(vcCrs[dir][ind[dir]]), 3)
              for ind in CartesianIndices(cells(rCrs))]
    @test first.(fstBlk) == vec(expPos)
    @test all(trp -> trp[2] ≈ Float64(prod(scale(rCrs))), fstBlk)
    # The second block reads off the fine region grid the same way
    vcFin = voxel_centers(rFin)
    expPos2 = [ntuple(dir -> Float64(vcFin[dir][ind[dir]]), 3)
               for ind in CartesianIndices(cells(rFin))]
    @test first.(sndBlk) == vec(expPos2)

    # cellvolumes matches the cells block by block
    celVol = cellvolumes(g)
    nCrs = prod(cells(rCrs))
    nFin = prod(cells(rFin))
    @test length(celVol) == nCrs + nFin
    @test all(celVol[1:nCrs] .≈ Float64(prod(scale(rCrs))))
    @test all(celVol[(nCrs + 1):end] .≈ Float64(prod(scale(rFin))))
    # The weighted sum over the cells is the domain volume
    @test sum(celVol) ≈ Float64(tiled_volume(g))
    # cellvolumes matches the volumes reported by coordinates
    @test celVol ≈ [trp[2] for trp in crd]

    # Region order matters for equality
    @test CompositeGrid([rCrs, rFin]) == CompositeGrid([rCrs, rFin])
    @test CompositeGrid([rCrs, rFin]) != CompositeGrid([rFin, rCrs])
    @test occursin("CompositeGrid (2 regions)", sprint(show, g))
end

@testset "CompositeGrid refine convenience form" begin
    r = mkgrid()
    box = (grdCtr0, (1 // 8, 1 // 8, 1 // 8))
    @test refine(r, box) == refine(CompositeGrid(r), box)
    @test refine(r, box; factor=(1, 1, 2)) ==
          refine(CompositeGrid(r), box; factor=(1, 1, 2))
    # A Region works as a box: only its outer bounds matter
    boxReg = Region((2, 2, 2), grdScl16, grdCtr0)
    @test refine(r, boxReg) == refine(r, (grdCtr0, (1 // 8, 1 // 8, 1 // 8)))
    # A finer box region with the same bounds gives the same tiling
    boxFin = Region((4, 4, 4), grdScl32, grdCtr0)
    @test refine(r, boxFin) == refine(r, boxReg)
end

@testset "CompositeGrid coordinates ordering" begin
    # coordinates must walk regions in order and, inside a region, column major
    g = refine(mkgrid(), (grdCtr0, (1 // 8, 1 // 8, 1 // 8)))
    expected = Tuple{NTuple{3,Float64},Float64,Int}[]
    for (idx, reg) in enumerate(regions(g))
        vc = voxel_centers(reg)
        vol = Float64(prod(scale(reg)))
        for ind in CartesianIndices(cells(reg))
            push!(expected, (ntuple(dir -> Float64(vc[dir][ind[dir]]), 3), vol, idx))
        end
    end
    @test collect(coordinates(g)) == expected
    @test length(expected) == length(g)
    # Total volume from the iterator equals the exact rational domain volume
    @test sum(trp -> trp[2], expected) ≈ Float64(tiled_volume(g))
end
