# Anti-aliasing strategies, exercised through `aa(shape, vc, vs, bg, strategy)`.
# `aa` returns the background when the voxel is outside the shape, the fill value
# when fully inside, and a blended value on the boundary.

@testset "AntiAliasing" begin
    sphere = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)

    @testset "NoAntiAliasing is a hard point test at the voxel center" begin
        @test aa(sphere, (0.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, NoAntiAliasing()) == 1.0
        @test aa(sphere, (5.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, NoAntiAliasing()) == 0.0
    end

    @testset "SuperResolutionAntiAliasing averages sub-samples" begin
        # Half-space z<=0. A voxel centered on the plane is half covered.
        h = FillableHalfSpace((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0)
        v = aa(h, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, SuperResolutionAntiAliasing(2))
        @test 0.0 < v < 1.0
        @test v ≈ 0.5
        # A voxel well inside is fully covered; well outside is empty.
        @test aa(h, (0.0, 0.0, -5.0), (1.0, 1.0, 1.0), 0.0, SuperResolutionAntiAliasing(4)) ≈ 1.0
        @test aa(h, (0.0, 0.0, 5.0), (1.0, 1.0, 1.0), 0.0, SuperResolutionAntiAliasing(4)) ≈ 0.0
        # Isotropic constructor with default n
        @test SuperResolutionAntiAliasing() isa SuperResolutionAntiAliasing
        @test SuperResolutionAntiAliasing(3).super_grid == (3, 3, 3)
        @test SuperResolutionAntiAliasing((2, 3, 4)).super_grid == (2, 3, 4)
    end

    @testset "SubpixelAntiAliasing uses the SDF for analytic coverage" begin
        h = FillableHalfSpace((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0)
        # On the boundary the coverage fraction is 1/2.
        @test aa(h, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, SubpixelAntiAliasing()) ≈ 0.5
        # Clearly inside -> fill; clearly outside -> background.
        @test aa(h, (0.0, 0.0, -5.0), (1.0, 1.0, 1.0), 0.0, SubpixelAntiAliasing()) ≈ 1.0
        @test aa(h, (0.0, 0.0, 5.0), (1.0, 1.0, 1.0), 0.0, SubpixelAntiAliasing()) == 0.0
    end

    @testset "GaussianAntiAliasing" begin
        # Kernel size must be odd.
        @test_throws AssertionError GaussianAntiAliasing((1.0, 1.0, 1.0), (2, 3, 3))
        big = FillableCube((0.0, 0.0, 0.0), 200.0, 1.0)  # huge box
        g = GaussianAntiAliasing(1.0, 3)
        # Stencil fully inside -> normalized weights sum to fill value.
        @test aa(big, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, g) ≈ 1.0
        # Stencil fully outside -> background.
        @test aa(big, (1000.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0, g) ≈ 0.0
    end

    @testset "AdaptiveAntiAliasing skips the stencil away from the boundary" begin
        cube = FillableCube((0.0, 0.0, 0.0), 2.0, 1.0)   # exact SDF
        a = AdaptiveAntiAliasing(SuperResolutionAntiAliasing(2))
        # Far inside -> fill, far outside -> background, via the single SDF check.
        @test aa(cube, (0.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, a) ≈ 1.0
        @test aa(cube, (5.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, a) ≈ 0.0
        # For a shape without an exact SDF it always defers to the inner strategy.
        sph = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)
        @test aa(sph, (0.0, 0.0, 0.0), (0.1, 0.1, 0.1), 0.0, a) ≈ 1.0
    end
end
