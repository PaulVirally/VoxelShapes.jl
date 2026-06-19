# Fill functions are callables taking a 3-tuple of local coordinates and
# returning a fill value. Tested directly as callables here; their use inside
# shapes is covered in the integration pass.

@testset "Fills" begin
    @testset "ConstantFill" begin
        f = ConstantFill(3.5)
        @test f((0.0, 0.0, 0.0)) == 3.5
        @test f((1.0, -2.0, 0.5)) == 3.5
        # works with non-scalar values too
        g = ConstantFill(SVector(1.0, 2.0, 3.0))
        @test g((0.1, 0.2, 0.3)) == SVector(1.0, 2.0, 3.0)
    end

    @testset "RadialGradient" begin
        # inner_value at r=0, outer_value at r=1 (r = norm(local), clamped to [0,1])
        f = RadialGradient(1.0, 0.0)
        @test f((0.0, 0.0, 0.0)) ≈ 1.0          # center
        @test f((1.0, 0.0, 0.0)) ≈ 0.0          # surface
        @test f((0.5, 0.0, 0.0)) ≈ 0.5          # halfway
        # r is clamped: beyond the surface stays at outer_value
        @test f((2.0, 0.0, 0.0)) ≈ 0.0
        # radius combines all three local coordinates
        @test f((0.0, 0.6, 0.8)) ≈ 0.0          # norm = 1
    end

    @testset "AxialGradient" begin
        # interpolates along a chosen local axis from -1 (v0) to +1 (v1)
        f = AxialGradient(3, 0.0, 1.0)
        @test f((0.0, 0.0, -1.0)) ≈ 0.0         # v0 end
        @test f((0.0, 0.0, 1.0)) ≈ 1.0          # v1 end
        @test f((0.0, 0.0, 0.0)) ≈ 0.5          # midpoint
        # clamped outside [-1, 1]
        @test f((0.0, 0.0, -5.0)) ≈ 0.0
        @test f((0.0, 0.0, 5.0)) ≈ 1.0
        # axis selection: interpolate along x instead
        fx = AxialGradient(1, 10.0, 20.0)
        @test fx((-1.0, 0.0, 0.0)) ≈ 10.0
        @test fx((1.0, 0.0, 0.0)) ≈ 20.0
    end
end
