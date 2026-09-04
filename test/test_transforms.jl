# Rotated wraps any shape and maps world points into the shape's local frame
# before the containment test. Rotation is isometric, so SDFs are preserved and
# has_exact_sdf / interpolation delegate to the inner shape.

@testset "Transforms (Rotated)" begin
    # An asymmetric box: long along x (half-length 3), thin along y and z.
    box = FillableCuboid((0.0, 0.0, 0.0), (6.0, 2.0, 2.0), 1.0)

    @testset "identity rotations leave the shape unchanged" begin
        # axis-angle with zero angle
        r0 = Rotated(box, (0.0, 0.0, 1.0), 0.0)
        @test ((2.5, 0.0, 0.0) in r0) == ((2.5, 0.0, 0.0) in box)
        @test ((0.0, 2.5, 0.0) in r0) == ((0.0, 2.5, 0.0) in box)
        # explicit identity matrix
        I3 = SMatrix{3,3,Float64,9}(1,0,0, 0,1,0, 0,0,1)
        rI = Rotated(box, I3)
        @test ((2.5, 0.0, 0.0) in rI)
        @test !((0.0, 2.5, 0.0) in rI)
        # Euler angles all zero
        rE = Rotated(box, (0.0, 0.0, 0.0))
        @test ((2.5, 0.0, 0.0) in rE)
    end

    @testset "90 degrees about z swaps the long axis from x to y" begin
        r = Rotated(box, (0.0, 0.0, π/2))   # Euler ZYX, only αz
        # Before rotation the long axis is x: (2.5,0,0) in, (0,2.5,0) out.
        # After a 90 degree turn about z, the long axis points along world y.
        @test (0.0, 2.5, 0.0) in r
        @test !((2.5, 0.0, 0.0) in r)
    end

    @testset "axis-angle about z, 90 degrees" begin
        r = Rotated(box, (0.0, 0.0, 1.0), π/2)
        @test (0.0, 2.5, 0.0) in r
        @test !((2.5, 0.0, 0.0) in r)
    end

    @testset "delegation to inner shape" begin
        r = Rotated(box, (0.0, 0.0, π/3))
        @test interpolation(r) === interpolation(box)
        @test has_exact_sdf(r) == has_exact_sdf(box)   # box is exact -> true
        # non-exact inner shape stays non-exact under rotation
        sph = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)
        @test has_exact_sdf(Rotated(sph, (0.0, 0.0, 0.3))) == false
    end

    @testset "rotation preserves distance (isometry)" begin
        r = Rotated(box, (0.0, 0.0, π/2))
        # The SDF of the rotated box at a world point equals the inner SDF at the
        # mapped point; magnitude of distance is preserved by rotation.
        # A point 1 unit beyond the (now-y-facing) long face:
        @test sdf(r, (0.0, 4.0, 0.0)) ≈ sdf(box, (4.0, 0.0, 0.0))
    end

    @testset "explicit pivot" begin
        # Rotate about a pivot away from the box center.
        I3 = SMatrix{3,3,Float64,9}(1,0,0, 0,1,0, 0,0,1)
        r = Rotated(box, I3, (10.0, 10.0, 10.0))
        # identity rotation about any pivot is still the identity
        @test (2.5, 0.0, 0.0) in r
    end

    @testset "bounding_box" begin
        @testset "identity rotation matches the inner box" begin
            I3 = SMatrix{3,3,Float64,9}(1,0,0, 0,1,0, 0,0,1)
            r = Rotated(box, I3)
            @test bounding_box(r) == bounding_box(box)
        end

        @testset "45 degrees about z widens x/y by sqrt(2)" begin
            # box has full side lengths (6,2,2) -> half-lengths (3,1,1) centered at the origin.
            r = Rotated(box, (0.0, 0.0, π/4))
            lo, hi = bounding_box(r)
            expected_xy = (3.0 + 1.0) / sqrt(2)   # widened half-extent along x and y
            @test lo[1] ≈ -expected_xy
            @test hi[1] ≈ expected_xy
            @test lo[2] ≈ -expected_xy
            @test hi[2] ≈ expected_xy
            @test lo[3] ≈ -1.0    # z untouched by a rotation about z
            @test hi[3] ≈ 1.0
        end

        @testset "conservativeness under random rotations" begin
            for _ in 1:200
                angles = (2π*rand(), 2π*rand(), 2π*rand())
                r = Rotated(box, angles)
                lo, hi = bounding_box(r)
                for _ in 1:20
                    p = ntuple(i -> 10.0 * (2*rand() - 1), 3)
                    if p in r
                        @test all(lo[i] <= p[i] <= hi[i] for i in 1:3)
                    end
                end
            end
        end

        @testset "infinite inner box stays all-infinite" begin
            # FillableHalfSpace has no center(), so build the rotation matrix
            # and pivot explicitly rather than going through the angle-based
            # constructors (which look up center(inner)).
            h = FillableHalfSpace((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0)
            c, s = cos(π/4), sin(π/4)
            Rz = SMatrix{3,3,Float64,9}(c, s, 0, -s, c, 0, 0, 0, 1)
            r = Rotated(h, Rz, (0.0, 0.0, 0.0))
            lo, hi = bounding_box(r)
            @test all(isinf, lo)
            @test all(isinf, hi)
        end
    end
end
