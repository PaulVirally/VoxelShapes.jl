# CompositeField: rasterizing a geometry into a composite grid, the flat layout,
# the per-region views, regrid, and refine driven by a shape's bounding box.

@testset "CompositeField" begin
    # An 8x8x8 volume at 1/16 refined by 2 inside a box of side 1/8 about the
    # origin gives 7 regions.
    refgrid() = refine(Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16)),
                       ((0 // 1, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8)))

    smallsphere() = FillableSphere((0.0, 0.0, 0.0), 0.05, 1.0)

    @testset "flat layout" begin
        grid = refgrid()
        field = rasterize(Geometry([smallsphere()], 0.0), grid)
        @test field isa CompositeField{Float64,Vector{Float64}}
        @test eltype(field) == Float64
        @test length(field) == length(grid)
        @test length(vec(field)) == length(grid)
        @test vec(field) === field.data
        @test regions(field) == regions(grid)
        @test nregions(field) == 7

        # regionview shapes match the region cell counts, and the views tile
        # the flat vector in region order.
        blocks = collect(eachregion(field))
        @test length(blocks) == nregions(grid)
        @test [size(b) for b in blocks] == [cells(r) for r in regions(grid)]
        @test vcat([vec(collect(b)) for b in blocks]...) == vec(field)
        for idx in 1:nregions(grid)
            @test size(regionview(field, idx)) == cells(regions(grid)[idx])
        end
    end

    @testset "each region block matches a standalone rasterize" begin
        geometry = Geometry([smallsphere()], 0.0)
        grid = refgrid()
        field = rasterize(geometry, grid)
        for (idx, reg) in enumerate(regions(grid))
            @test collect(regionview(field, idx)) == rasterize(geometry, reg)
        end
    end

    @testset "sample agrees with the field at every voxel center" begin
        geometry = Geometry([smallsphere()], 0.0)
        grid = refgrid()
        field = rasterize(geometry, grid)
        for (flat, (pos, _, regidx)) in enumerate(coordinates(grid))
            vs = Float64.(scale(regions(grid)[regidx]))
            @test vec(field)[flat] == sample(geometry, pos, vs)
        end
    end

    @testset "regrid onto the finest scale" begin
        geometry = Geometry([smallsphere()], 0.0)
        grid = refgrid()
        field = rasterize(geometry, grid)
        uniform = regrid(field)
        lwr, upr = bounding_box(grid)
        fine = Region((16, 16, 16), finest(grid),
                      ntuple(dir -> (lwr[dir] + upr[dir]) // 2, 3))
        @test size(uniform) == (16, 16, 16)
        @test size(uniform) == cells(fine)

        # The sphere fits inside the refined core, so every other region is
        # empty and its coarse cells replicate to zero. The uniform array
        # equals a plain rasterize at the finest scale.
        dense = rasterize(geometry, fine)
        @test uniform == dense

        # The core block is exactly the matching slice of the dense array.
        core = regionview(field, 1)
        @test size(core) == (8, 8, 8)
        @test collect(core) == dense[5:12, 5:12, 5:12]

        # regrid defaults to the finest scale, and an explicit scale agrees.
        @test regrid(field, finest(grid)) == uniform
    end

    @testset "regrid replicates coarse voxels" begin
        # A shape that covers a whole coarse region shows up as a solid block of
        # finest voxels after regridding.
        grid = refgrid()
        cube = FillableCuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 3.0)
        field = rasterize(Geometry([cube], 0.0), grid)
        uniform = regrid(field)
        @test all(==(3.0), uniform)
        @test sum(uniform) == 3.0 * 16^3
    end

    @testset "regrid rejects a scale that does not divide the grid" begin
        field = rasterize(Geometry([smallsphere()], 0.0), refgrid())
        # 1/48 does not divide 1/32, the scale of the refined core.
        @test_throws ArgumentError regrid(field, (1 // 48, 1 // 48, 1 // 48))
        # 1/16 divides neither the core scale nor, going the other way, the box.
        @test_throws ArgumentError regrid(field, (1 // 16, 1 // 16, 1 // 16))
    end

    @testset "coarser scales that do divide are allowed" begin
        # A one-region grid at 1/8 regridded onto 1/8 is the region itself.
        grid = CompositeGrid(Region((2, 2, 2), (1 // 8, 1 // 8, 1 // 8)))
        geometry = Geometry([FillableSphere((0.0, 0.0, 0.0), 0.2, 1.0)], 0.0)
        field = rasterize(geometry, grid)
        @test regrid(field) == rasterize(geometry, regions(grid)[1])
    end

    @testset "equality and show" begin
        geometry = Geometry([smallsphere()], 0.0)
        grid = refgrid()
        a = rasterize(geometry, grid)
        b = rasterize(geometry, grid)
        @test a == b
        c = rasterize(Geometry([smallsphere()], 1.0), grid)
        @test a != c
        # Different grid, same length: not equal.
        other = refine(Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16)),
                       ((0 // 1, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8));
                       factor=(2, 2, 2))
        @test rasterize(geometry, other) == a   # same grid built the same way
        str = sprint(show, a)
        @test occursin("CompositeField", str)
        @test occursin("7 regions", str)
        @test occursin(string(length(grid)), str)
    end

    @testset "CompositeField constructor checks the length" begin
        grid = refgrid()
        @test_throws DimensionMismatch CompositeField(zeros(3), grid)
        f = CompositeField(collect(1.0:length(grid)), grid)
        @test f.offsets[1] == 0
        @test f.offsets[end] == length(grid)
        @test f.offsets == cumsum([0; [length(r) for r in regions(grid)]])
        # The block boundaries agree with the region views.
        @test regionview(f, 1)[1] == 1.0
        @test regionview(f, 2)[1] == length(regions(grid)[1]) + 1.0
    end

    @testset "rasterize! refills an existing field" begin
        grid = refgrid()
        field = CompositeField(fill(9.0, length(grid)), grid)
        geometry = Geometry([smallsphere()], 0.0)
        @test rasterize!(field, geometry) === field
        @test field == rasterize(geometry, grid)
    end

    @testset "refine around a shape puts its bounding box in the core" begin
        base = Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16))
        sphere = FillableSphere((0.0, 0.0, 0.0), 0.05, 1.0)
        grid = refine(base, sphere)
        @test grid == refgrid()

        # Every point of the bounding box lies inside the refined core.
        core = regions(grid)[1]
        lwr, upr = bounding_box(sphere)
        @test all(lower_corner(core) .<= lwr)
        @test all(upr .<= upper_corner(core))
        @test scale(core) == (1 // 32, 1 // 32, 1 // 32)

        # Padding grows the box, so the core grows with it.
        padded = refine(base, sphere; padding=0.1)
        pcore = regions(padded)[1]
        @test all(lower_corner(pcore) .<= lower_corner(core))
        @test all(upper_corner(pcore) .>= upper_corner(core))
        @test all(lower_corner(pcore) .<= lwr .- 0.1)
        @test all(upr .+ 0.1 .<= upper_corner(pcore))
        # A per-axis padding is allowed too.
        @test refine(base, sphere; padding=(0.1, 0.1, 0.1)) == padded
    end

    @testset "refine around a shape that misses the grid" begin
        base = CompositeGrid(Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16)))
        far = FillableSphere((100.0, 0.0, 0.0), 0.05, 1.0)
        @test refine(base, far) === base
        # A shape whose box only touches the grid face has empty interior there.
        touching = FillableCuboid((-0.5, 0.0, 0.0), (0.5, 0.1, 0.1), 1.0)
        @test bounding_box(touching)[2][1] == -0.25
        @test refine(base, touching) === base
    end

    @testset "refine around a shape with no bounding box refines everything" begin
        # A half space with a diagonal normal is unbounded on every axis. The
        # clamped box is the whole grid, so every region is carved.
        base = Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16))
        half = FillableHalfSpace((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 1.0)
        grid = refine(base, half)
        @test nregions(grid) == 1
        @test cells(regions(grid)[1]) == (16, 16, 16)
        @test scale(regions(grid)[1]) == (1 // 32, 1 // 32, 1 // 32)
    end

    @testset "refine around a shape rejects bad arguments" begin
        base = Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16))
        sphere = FillableSphere((0.0, 0.0, 0.0), 0.05, 1.0)
        @test_throws ArgumentError refine(base, sphere; snap=:inward)
        @test_throws ArgumentError refine(base, sphere; factor=0)
        @test_throws ArgumentError refine(base, sphere; factor=1.5)
    end
end
