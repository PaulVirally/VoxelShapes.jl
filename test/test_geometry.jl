# Geometry: what to draw, and single region rasterization onto a Region. The
# constructors, eltype, add_shape, sample, rasterize and rasterize!.

@testset "Geometry" begin
    # A 3x3x3 grid of unit voxels with its lower corner at the origin has its
    # center at (3/2, 3/2, 3/2).
    unitgrid() = Region((3, 3, 3), (1 // 1, 1 // 1, 1 // 1), (3 // 2, 3 // 2, 3 // 2))

    @testset "construction and interface" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        sv = Geometry([sphere], 0.0, NoAntiAliasing())
        @test eltype(sv) == Float64
        @test sv.shapes isa Tuple
        @test length(sv.shapes) == 1
        st = Geometry((sphere,), 0.0, NoAntiAliasing())
        @test st.shapes == (sphere,)
        # The anti-aliasing strategy defaults to NoAntiAliasing.
        @test Geometry((sphere,), 0.0).aa == NoAntiAliasing()
        @test Geometry([sphere], 0.0).aa == NoAntiAliasing()
        # A geometry of isbits shapes is itself isbits. The GPU kernel needs
        # this.
        @test isbits(st)
    end

    @testset "eltype follows the background value" begin
        sphere = FillableSphere((0.0, 0.0, 0.0), 1.0, 1)
        @test eltype(Geometry([sphere], 0)) == Int
        @test eltype(Geometry([sphere], 0.0 + 0.0im)) == ComplexF64
    end

    @testset "add_shape appends and is non-mutating" begin
        a = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)
        b = FillableSphere((2.0, 0.0, 0.0), 1.0, 2.0)
        s = Geometry([a], 0.0)
        s2 = add_shape(s, b)
        @test length(s.shapes) == 1       # original untouched
        @test length(s2.shapes) == 2
        @test s2 isa Geometry
    end

    @testset "rasterize over a Region" begin
        # A small sphere at (1.5,1.5,1.5) only covers the center voxel (2,2,2).
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        arr = rasterize(Geometry([sphere], 0.0), unitgrid())
        @test arr isa Array{Float64,3}
        @test size(arr) == (3, 3, 3)
        @test arr[2, 2, 2] == 1.0
        @test sum(arr) == 1.0             # exactly one voxel filled
    end

    @testset "rasterize! into a preallocated array matches rasterize" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        geometry = Geometry([sphere], 0.0)
        arr = zeros(Float64, 3, 3, 3)
        @test rasterize!(arr, geometry, unitgrid()) === arr
        @test arr == rasterize(geometry, unitgrid())
        @test_throws DimensionMismatch rasterize!(zeros(2, 2, 2), geometry, unitgrid())
    end

    @testset "empty geometry is all background" begin
        arr = rasterize(Geometry((), 7.0), Region((2, 2, 2), (1 // 1, 1 // 1, 1 // 1)))
        @test all(==(7.0), arr)
    end

    @testset "shape ordering: first match claims the voxel" begin
        # Two cubes both covering the single voxel of a 1x1x1 grid.
        c1 = FillableCube((0.5, 0.5, 0.5), 10.0, 1.0)
        c2 = FillableCube((0.5, 0.5, 0.5), 10.0, 2.0)
        reg = Region((1, 1, 1), (1 // 1, 1 // 1, 1 // 1), (1 // 2, 1 // 2, 1 // 2))
        # first listed wins
        @test rasterize(Geometry([c1, c2], 0.0), reg)[1, 1, 1] == 1.0
        @test rasterize(Geometry([c2, c1], 0.0), reg)[1, 1, 1] == 2.0
    end

    @testset "the region center places the grid in space" begin
        # A region centered on (11,11,11) with 2 unit voxels per axis has its
        # first voxel center at (10.5,10.5,10.5).
        sphere = FillableSphere((10.5, 10.5, 10.5), 0.3, 1.0)
        reg = Region((2, 2, 2), (1 // 1, 1 // 1, 1 // 1), (11 // 1, 11 // 1, 11 // 1))
        arr = rasterize(Geometry([sphere], 0.0), reg)
        @test arr[1, 1, 1] == 1.0
        @test sum(arr) == 1.0
    end

    @testset "voxel centers follow the region" begin
        # Rasterizing a geometry must agree with sampling it at voxel_centers.
        sphere = FillableSphere((0.0, 0.0, 0.0), 0.9, 1.0)
        geometry = Geometry([sphere], 0.0, SuperResolutionAntiAliasing(2))
        reg = Region((4, 4, 4), (1 // 2, 1 // 2, 1 // 2))
        arr = rasterize(geometry, reg)
        vcs = voxel_centers(reg)
        vs = Float64.(scale(reg))
        for idx in CartesianIndices(arr)
            pos = ntuple(dir -> Float64(vcs[dir][idx[dir]]), 3)
            @test arr[idx] == sample(geometry, pos, vs)
        end
    end

    @testset "precision controls the kernel's float type" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 1.0f0)
        arr = rasterize(Geometry([sphere], 0.0f0), unitgrid(); precision=Float32)
        @test arr isa Array{Float32,3}
        @test arr[2, 2, 2] == 1.0f0
        @test sum(arr) == 1.0f0
    end
end
