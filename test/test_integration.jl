# Pass 2: cross-cutting behavior that the per-feature tests above don't reach.
# Targets the seams between components.

@testset "Integration" begin
    # A 3x3x3 grid of unit voxels with its lower corner at the origin.
    unitgrid() = Region((3, 3, 3), (1 // 1, 1 // 1, 1 // 1), (3 // 2, 3 // 2, 3 // 2))

    @testset "interpolation flows through anti-aliasing" begin
        # On a boundary voxel, the shape's interpolation decides how
        # sub-samples combine. Max claims the voxel if any sub-sample is
        # inside, min only if all of them are.
        plane = (0.0, 0.0, 1.0)
        h_max = FillableHalfSpace((0.0, 0.0, 0.0), plane, 1.0; interpolation=MaxInterpolation())
        h_min = FillableHalfSpace((0.0, 0.0, 0.0), plane, 1.0; interpolation=MinInterpolation())
        sr = SuperResolutionAntiAliasing(2)
        # Voxel straddling the plane: some sub-samples in, some out.
        @test aa(h_max, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, sr) ≈ 1.0
        @test aa(h_min, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, sr) ≈ 0.0
    end

    @testset "AdaptiveAntiAliasing is safe on shapes without an SDF" begin
        # A cone has no exact SDF. Adaptive must not call sdf on it and
        # instead defers to the inner strategy.
        cone = FillableCone((0.0, 0.0, 0.0), 1.0, 0.0, 1.0, 1.0)
        a = AdaptiveAntiAliasing(SuperResolutionAntiAliasing(2))
        @test aa(cone, (0.0, 0.0, -0.9), (0.05, 0.05, 0.05), 0.0, a) ≈ 1.0

        # SubpixelAntiAliasing requires an SDF, so it errors on a cone.
        @test_throws MethodError aa(cone, (0.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, SubpixelAntiAliasing())
    end

    @testset "rotating a shape that has no center needs an explicit pivot" begin
        cap = FillableCapsule((0.0, 0.0, 0.0), (0.0, 0.0, 2.0), 0.5, 1.0)
        # Capsule defines no center(), so the pivot-inferring constructors fail.
        @test_throws MethodError Rotated(cap, (0.0, 0.0, 0.3))
        # Supplying an explicit pivot works.
        I3 = SMatrix{3,3,Float64,9}(1,0,0, 0,1,0, 0,0,1)
        r = Rotated(cap, I3, (0.0, 0.0, 1.0))
        @test (0.0, 0.0, 1.0) in r
    end

    @testset "the background value is a transparency sentinel" begin
        # A shape that evaluates to the background value counts as
        # transparent, so evaluation falls through to the next shape.
        c1 = FillableCube((0.5, 0.5, 0.5), 10.0, 0.0)   # fill equals background
        c2 = FillableCube((0.5, 0.5, 0.5), 10.0, 5.0)
        reg = Region((1, 1, 1), (1 // 1, 1 // 1, 1 // 1), (1 // 2, 1 // 2, 1 // 2))
        # c1 is transparent, so c2 shows through.
        @test rasterize(Geometry([c1, c2], 0.0), reg)[1, 1, 1] == 5.0

        # With no lower layer to show through, the voxel stays the background.
        @test rasterize(Geometry([c1], 0.0), reg)[1, 1, 1] == 0.0
    end

    @testset "gradient fill rasterized through a Geometry" begin
        # The convenience constructors wrap fill_val in a closure. A gradient
        # must go through the struct's inner constructor as fill_function.
        g = RadialGradient(1.0, 0.0)
        shape = FillableEllipsoid{Float64, typeof(g), LinearInterpolation}(
            SVector(1.5, 1.5, 1.5), SVector(0.6, 0.6, 0.6), g, LinearInterpolation())
        arr = rasterize(Geometry([shape], 0.0), unitgrid())
        # The center voxel is at the sphere center: local coords (0,0,0).
        @test arr[2, 2, 2] ≈ 1.0
    end

    @testset "CSG rasterized through a Geometry" begin
        a = FillableCube((1.5, 1.5, 1.5), 1.2, 1.0)
        b = FillableCube((1.5, 1.5, 1.5), 1.2, 2.0)
        u = csg_union(a, b)
        arr = rasterize(Geometry([u], 0.0), unitgrid())
        # Center voxel is inside a -> union fill delegates to a.
        @test arr[2, 2, 2] == 1.0
    end

    @testset "Rotated shape rasterized through a Geometry" begin
        # A box long along x, thin in y/z. Rotating 90 deg about z reorients it
        # along y, so the filled voxels move from a row to a column.
        box = FillableCuboid((2.5, 2.5, 2.5), (4.0, 0.8, 0.8), 1.0)
        reg = Region((5, 5, 5), (1 // 1, 1 // 1, 1 // 1), (5 // 2, 5 // 2, 5 // 2))
        rot = Rotated(box, (0.0, 0.0, π/2))
        au = rasterize(Geometry([box], 0.0), reg)
        ar = rasterize(Geometry([rot], 0.0), reg)

        @test au[1, 3, 3] == 1.0   # unrotated: filled along x
        @test au[3, 1, 3] == 0.0
        @test ar[3, 1, 3] == 1.0   # rotated: filled along y
        @test ar[1, 3, 3] == 0.0
        @test sum(au) == sum(ar)   # rotation preserves the filled count
    end

    @testset "vector-valued fill through a Geometry with interpolation" begin
        bg = SVector(0.0, 0.0, 0.0)
        s = FillableSphere((1.5, 1.5, 1.5), 0.6, SVector(1.0, 2.0, 3.0))
        arr = rasterize(Geometry([s], bg), unitgrid())
        @test eltype(arr) == SVector{3, Float64}
        @test arr[2, 2, 2] == SVector(1.0, 2.0, 3.0)
        @test arr[1, 1, 1] == bg
    end

    @testset "anti-aliased geometry produces partial-coverage values" begin
        # A sphere boundary with super-resolution yields voxels strictly
        # between background and fill somewhere on the surface.
        s = FillableSphere((1.5, 1.5, 1.5), 1.0, 1.0)
        arr = rasterize(Geometry([s], 0.0, SuperResolutionAntiAliasing(4)), unitgrid())
        @test any(v -> 0.0 < v < 1.0, arr)     # boundary blending occurred
        @test all(v -> 0.0 <= v <= 1.0, arr)   # stays within [bg, fill]
    end
end
