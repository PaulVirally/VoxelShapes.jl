# World: the immutable scene descriptor, its constructors, the array-like
# interface (size/eltype/similar), add_shape, and CPU rasterization via
# Array(world) / fill!.

@testset "World" begin
    @testset "construction and interface" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        # Vector of shapes (convenience constructor), default origin.
        w = World((3, 3, 3), (1.0, 1.0, 1.0), [sphere], 0.0, NoAntiAliasing())
        @test size(w) == (3, 3, 3)
        @test eltype(w) == Float64
        # Tuple of shapes constructor.
        wt = World((3, 3, 3), (1.0, 1.0, 1.0), (sphere,), 0.0, NoAntiAliasing())
        @test size(wt) == (3, 3, 3)
        # similar gives a zeroed array of the right size and element type.
        s = similar(w)
        @test size(s) == (3, 3, 3)
        @test eltype(s) == Float64
        @test all(iszero, s)
    end

    @testset "eltype follows the background value" begin
        sphere = FillableSphere((0.0, 0.0, 0.0), 1.0, 1)
        wi = World((2, 2, 2), (1.0, 1.0, 1.0), [sphere], 0, NoAntiAliasing())
        @test eltype(wi) == Int
        wc = World((2, 2, 2), (1.0, 1.0, 1.0), [sphere], 0.0 + 0.0im, NoAntiAliasing())
        @test eltype(wc) == ComplexF64
    end

    @testset "add_shape appends and is non-mutating" begin
        a = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)
        b = FillableSphere((2.0, 0.0, 0.0), 1.0, 2.0)
        w = World((3, 3, 3), (1.0, 1.0, 1.0), [a], 0.0, NoAntiAliasing())
        w2 = add_shape(w, b)
        @test length(w.shapes) == 1       # original untouched
        @test length(w2.shapes) == 2
        @test w2 isa World
    end

    @testset "Array(world) rasterizes on the CPU" begin
        # 3x3x3 grid, unit voxels, origin 0. Voxel (i,j,k) center = (i-0.5, j-0.5, k-0.5).
        # A small sphere at (1.5,1.5,1.5) only covers the center voxel (2,2,2).
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        w = World((3, 3, 3), (1.0, 1.0, 1.0), [sphere], 0.0, NoAntiAliasing())
        arr = Array(w)
        @test arr isa Array{Float64,3}
        @test size(arr) == (3, 3, 3)
        @test arr[2, 2, 2] == 1.0
        @test sum(arr) == 1.0             # exactly one voxel filled
    end

    @testset "fill! into a preallocated array matches Array(world)" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0)
        w = World((3, 3, 3), (1.0, 1.0, 1.0), [sphere], 0.0, NoAntiAliasing())
        arr = zeros(Float64, 3, 3, 3)
        fill!(arr, w)
        @test arr == Array(w)
    end

    @testset "empty world is all background" begin
        w = World((2, 2, 2), (1.0, 1.0, 1.0), (), 7.0, NoAntiAliasing())
        @test all(==(7.0), Array(w))
    end

    @testset "shape ordering: first match claims the voxel" begin
        # Two cubes both covering the center voxel of a 1x1x1 grid.
        c1 = FillableCube((0.5, 0.5, 0.5), 10.0, 1.0)
        c2 = FillableCube((0.5, 0.5, 0.5), 10.0, 2.0)
        w12 = World((1, 1, 1), (1.0, 1.0, 1.0), [c1, c2], 0.0, NoAntiAliasing())
        w21 = World((1, 1, 1), (1.0, 1.0, 1.0), [c2, c1], 0.0, NoAntiAliasing())
        @test Array(w12)[1, 1, 1] == 1.0   # first listed wins
        @test Array(w21)[1, 1, 1] == 2.0
    end

    @testset "origin offsets the grid in world space" begin
        # With origin (10,10,10) the center of voxel (1,1,1) is (10.5,10.5,10.5).
        sphere = FillableSphere((10.5, 10.5, 10.5), 0.3, 1.0)
        w = World((2, 2, 2), (1.0, 1.0, 1.0), [sphere], 0.0, NoAntiAliasing(),
                  (10.0, 10.0, 10.0))
        arr = Array(w)
        @test arr[1, 1, 1] == 1.0
        @test sum(arr) == 1.0
    end
end
