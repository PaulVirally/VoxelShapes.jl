# Per-primitive behavior: containment (`in`), fill value (`fill`), accessors,
# signed distance (`sdf`) and the `has_exact_sdf` flag. Geometry conventions are
# taken from the docstrings / README, not the implementation.

const VS = 1.0  # nominal voxel size used where fill needs a voxel_size argument

# A minimal shape with no bounding_box method, to exercise the Types default.
struct _DummyShape <: AbstractFillableShape end

@testset "Shapes" begin

    @testset "FillableSphere / FillableEllipsoid" begin
        s = FillableSphere((0.0, 0.0, 0.0), 1.0, 5.0)
        @test s isa AbstractFillableShape
        @test center(s) == SVector(0.0, 0.0, 0.0)
        @test radii(s) == SVector(1.0, 1.0, 1.0)
        @test interpolation(s) isa LinearInterpolation

        # containment
        @test (0.0, 0.0, 0.0) in s
        @test (0.99, 0.0, 0.0) in s
        @test (1.0, 0.0, 0.0) in s          # surface counts as inside
        @test !((1.01, 0.0, 0.0) in s)
        @test !((1.0, 1.0, 1.0) in s)

        # fill returns the constant value anywhere
        @test fill(s, (0.0, 0.0, 0.0), (VS, VS, VS)) == 5.0
        @test fill(s, (0.5, 0.0, 0.0), (VS, VS, VS)) == 5.0

        # approximate SDF: sign is correct (negative inside, positive outside)
        @test sdf(s, (0.5, 0.0, 0.0)) < 0
        @test sdf(s, (2.0, 0.0, 0.0)) > 0
        @test has_exact_sdf(s) == false

        # ellipsoid with independent radii
        e = FillableEllipsoid((1.0, 2.0, 3.0), (2.0, 1.0, 0.5), 1.0)
        @test radii(e) == SVector(2.0, 1.0, 0.5)
        @test center(e) == SVector(1.0, 2.0, 3.0)
        @test (1.0, 2.0, 3.0) in e               # center
        @test (3.0, 2.0, 3.0) in e               # surface tip along x (rx=2)
        @test !((1.0, 4.0, 3.0) in e)            # 2 beyond ry=1
        @test (1.0, 2.0, 3.45) in e              # within rz=0.5

        # custom interpolation keyword
        s2 = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0; interpolation=MaxInterpolation())
        @test interpolation(s2) isa MaxInterpolation
    end

    @testset "FillableCuboid / FillableCube" begin
        c = FillableCube((0.0, 0.0, 0.0), 2.0, 9.0)   # side length 2, half-length 1
        @test half_lengths(c) == SVector(1.0, 1.0, 1.0)
        @test lengths(c) == SVector(2.0, 2.0, 2.0)
        @test center(c) == SVector(0.0, 0.0, 0.0)

        @test (0.0, 0.0, 0.0) in c
        @test (1.0, 1.0, 1.0) in c               # corner is inside (inclusive)
        @test !((1.01, 0.0, 0.0) in c)
        @test fill(c, (0.3, 0.2, 0.1), (VS, VS, VS)) == 9.0

        # cuboid with independent full side lengths
        b = FillableCuboid((0.0, 0.0, 0.0), (4.0, 2.0, 1.0), 1.0)
        @test lengths(b) == SVector(4.0, 2.0, 1.0)
        @test half_lengths(b) == SVector(2.0, 1.0, 0.5)
        @test (2.0, 0.0, 0.0) in b
        @test !((2.01, 0.0, 0.0) in b)
        @test (0.0, 1.0, 0.5) in b

        # exact SDF: signed Euclidean distance to the box surface
        @test has_exact_sdf(c) == true
        @test sdf(c, (2.0, 0.0, 0.0)) ≈ 1.0      # 1 unit outside the +x face
        @test sdf(c, (0.0, 0.0, 0.0)) ≈ -1.0     # 1 unit inside, to nearest face
        @test sdf(c, (1.0, 0.0, 0.0)) ≈ 0.0      # on the surface
    end

    @testset "FillableCylinder" begin
        cyl = FillableCylinder((0.0, 0.0, 0.0), 1.0, 2.0, 3.0)  # axis=3, fill 3
        @test radius(cyl) == 1.0
        @test half_height(cyl) == 2.0
        @test center(cyl) == SVector(0.0, 0.0, 0.0)

        @test (0.0, 0.0, 0.0) in cyl
        @test (1.0, 0.0, 0.0) in cyl             # on the radial surface
        @test (0.0, 0.0, 2.0) in cyl             # on the end cap
        @test !((1.5, 0.0, 0.0) in cyl)          # outside radius
        @test !((0.0, 0.0, 2.5) in cyl)          # beyond half_height
        @test fill(cyl, (0.0, 0.0, 0.0), (VS, VS, VS)) == 3.0

        @test has_exact_sdf(cyl) == true
        @test sdf(cyl, (2.0, 0.0, 0.0)) ≈ 1.0    # 1 unit outside the curved face
        @test sdf(cyl, (0.0, 0.0, 3.0)) ≈ 1.0    # 1 unit beyond end cap
        @test sdf(cyl, (0.0, 0.0, 0.0)) < 0      # interior

        # axis selection: same cylinder oriented along x
        cylx = FillableCylinder((0.0, 0.0, 0.0), 1.0, 2.0, 1.0; axis=1)
        @test (2.0, 0.0, 0.0) in cylx            # axial extent now along x
        @test !((0.0, 0.0, 2.0) in cylx)         # not along z anymore
    end

    @testset "FillableHalfSpace" begin
        # inside = the side the normal points away from: dot(p - point, n) <= 0
        h = FillableHalfSpace((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 2.0)
        @test (0.0, 0.0, -1.0) in h
        @test (0.0, 0.0, 0.0) in h               # on the plane
        @test !((0.0, 0.0, 1.0) in h)
        @test fill(h, (0.0, 0.0, -1.0), (VS, VS, VS)) == 2.0

        @test has_exact_sdf(h) == true
        @test sdf(h, (0.0, 0.0, 1.0)) ≈ 1.0
        @test sdf(h, (0.0, 0.0, -1.0)) ≈ -1.0

        # the normal is normalized on construction
        hn = FillableHalfSpace((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), 1.0)
        @test sdf(hn, (0.0, 0.0, 2.0)) ≈ 2.0     # distance, not 2*5
    end

    @testset "FillableSlab" begin
        # half_thickness 1
        s = FillableSlab((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0, 7.0)
        @test (0.0, 0.0, 0.0) in s
        @test (0.0, 0.0, 1.0) in s               # on the slab face
        @test !((0.0, 0.0, 1.5) in s)
        @test (0.0, 0.0, -0.5) in s
        @test fill(s, (0.0, 0.0, 0.0), (VS, VS, VS)) == 7.0

        @test has_exact_sdf(s) == true
        @test sdf(s, (0.0, 0.0, 2.0)) ≈ 1.0      # 1 unit outside the slab
        # center, 1 unit from each face
        @test sdf(s, (0.0, 0.0, 0.0)) ≈ -1.0

        # normal normalization
        sn = FillableSlab((0.0, 0.0, 0.0), (0.0, 0.0, 3.0), 1.0, 1.0)
        @test sdf(sn, (0.0, 0.0, 3.0)) ≈ 2.0
    end

    @testset "FillableCone" begin
        # true cone: base_radius 1 at -half_height, top_radius 0 at +half_height
        cone = FillableCone((0.0, 0.0, 0.0), 1.0, 0.0, 1.0, 4.0)  # axis=3
        @test center(cone) == SVector(0.0, 0.0, 0.0)
        @test interpolation(cone) isa LinearInterpolation

        @test (0.0, 0.0, -1.0) in cone           # base center
        @test (0.9, 0.0, -1.0) in cone           # within base radius
        @test (0.0, 0.0, 1.0) in cone            # tip
        @test !((0.1, 0.0, 1.0) in cone)         # radius is 0 at the tip
        @test (0.5, 0.0, 0.0) in cone            # radius is 0.5 at mid-height
        @test !((0.6, 0.0, 0.0) in cone)
        @test !((0.0, 0.0, 1.5) in cone)         # beyond half_height
        @test fill(cone, (0.0, 0.0, 0.0), (VS, VS, VS)) == 4.0

        # frustum: nonzero unequal radii
        frust = FillableCone((0.0, 0.0, 0.0), 1.0, 0.5, 1.0, 1.0)
        @test (0.75, 0.0, 0.0) in frust          # mid radius is 0.75
        @test !((0.9, 0.0, 0.0) in frust)

        # no closed-form SDF
        @test has_exact_sdf(cone) == false
    end

    @testset "FillableTorus" begin
        t = FillableTorus((0.0, 0.0, 0.0), 2.0, 0.5, 6.0)  # R=2, r=0.5, axis=3
        @test major_radius(t) == 2.0
        @test minor_radius(t) == 0.5
        @test center(t) == SVector(0.0, 0.0, 0.0)

        @test (2.0, 0.0, 0.0) in t               # tube center
        @test (2.5, 0.0, 0.0) in t               # outer edge of tube
        @test (1.5, 0.0, 0.0) in t               # inner edge of tube
        @test !((3.0, 0.0, 0.0) in t)            # past the tube
        @test !((0.0, 0.0, 0.0) in t)            # the hole in the middle
        @test fill(t, (2.0, 0.0, 0.0), (VS, VS, VS)) == 6.0

        @test has_exact_sdf(t) == true
        @test sdf(t, (3.0, 0.0, 0.0)) ≈ 0.5
        @test sdf(t, (2.0, 0.0, 0.0)) ≈ -0.5
        @test sdf(t, (2.5, 0.0, 0.0)) ≈ 0.0
    end

    @testset "FillableCapsule" begin
        cap = FillableCapsule((0.0, 0.0, 0.0), (0.0, 0.0, 2.0), 0.5, 8.0)
        @test (0.0, 0.0, 1.0) in cap             # on the segment axis
        @test (0.4, 0.0, 1.0) in cap             # within tube radius
        @test !((0.6, 0.0, 1.0) in cap)          # outside tube radius
        @test (0.0, 0.0, 2.4) in cap             # inside the rounded end cap
        @test !((0.0, 0.0, 2.6) in cap)          # beyond the cap
        @test fill(cap, (0.0, 0.0, 1.0), (VS, VS, VS)) == 8.0

        @test has_exact_sdf(cap) == true
        @test sdf(cap, (1.0, 0.0, 1.0)) ≈ 0.5    # distance 1 to axis, - radius
        @test sdf(cap, (0.0, 0.0, 1.0)) ≈ -0.5
    end

    @testset "bounding_box" begin
        # Helper: draw random points in a cube and check that every point
        # inside the shape also lies inside the reported bounding box.
        function check_conservative(shape, sample_center, sample_halfwidth; n=2000)
            lo, hi = bounding_box(shape)
            for _ in 1:n
                p = ntuple(i -> sample_center[i] + sample_halfwidth[i] * (2*rand() - 1), 3)
                if p in shape
                    @test all(lo[i] <= p[i] for i in 1:3)
                    @test all(p[i] <= hi[i] for i in 1:3)
                end
            end
        end

        @testset "default method (unknown shape)" begin
            lo, hi = bounding_box(_DummyShape())
            @test lo == (-Inf, -Inf, -Inf)
            @test hi == (Inf, Inf, Inf)
        end

        @testset "Ellipsoid/Sphere" begin
            e = FillableEllipsoid((1.0, -2.0, 3.0), (2.0, 1.0, 0.5), 1.0)
            lo, hi = bounding_box(e)
            @test lo == (-1.0, -3.0, 2.5)
            @test hi == (3.0, -1.0, 3.5)
            check_conservative(e, (1.0, -2.0, 3.0), (3.0, 2.0, 1.5))

            s = FillableSphere((0.0, 0.0, 0.0), 2.0, 1.0)
            check_conservative(s, (0.0, 0.0, 0.0), (3.0, 3.0, 3.0))
        end

        @testset "Cuboid" begin
            c = FillableCuboid((1.0, 2.0, 3.0), (4.0, 2.0, 1.0), 1.0)
            lo, hi = bounding_box(c)
            @test lo == (-1.0, 1.0, 2.5)
            @test hi == (3.0, 3.0, 3.5)
            check_conservative(c, (1.0, 2.0, 3.0), (3.0, 2.5, 1.5))
        end

        @testset "Cylinder" begin
            cyl = FillableCylinder((0.0, 0.0, 0.0), 1.0, 2.0, 1.0; axis=3)
            lo, hi = bounding_box(cyl)
            @test lo == (-1.0, -1.0, -2.0)
            @test hi == (1.0, 1.0, 2.0)
            check_conservative(cyl, (0.0, 0.0, 0.0), (1.5, 1.5, 2.5))

            cylx = FillableCylinder((1.0, 0.0, 0.0), 0.5, 3.0, 1.0; axis=1)
            lo, hi = bounding_box(cylx)
            @test lo == (-2.0, -0.5, -0.5)
            @test hi == (4.0, 0.5, 0.5)
            check_conservative(cylx, (1.0, 0.0, 0.0), (3.5, 1.0, 1.0))
        end

        @testset "Cone" begin
            cone = FillableCone((0.0, 0.0, 0.0), 1.0, 0.5, 2.0, 1.0; axis=3)
            lo, hi = bounding_box(cone)
            @test lo == (-1.0, -1.0, -2.0)
            @test hi == (1.0, 1.0, 2.0)
            check_conservative(cone, (0.0, 0.0, 0.0), (1.5, 1.5, 2.5))
        end

        @testset "Torus" begin
            t = FillableTorus((0.0, 0.0, 0.0), 2.0, 0.5, 1.0; axis=3)
            lo, hi = bounding_box(t)
            @test lo == (-2.5, -2.5, -0.5)
            @test hi == (2.5, 2.5, 0.5)
            check_conservative(t, (0.0, 0.0, 0.0), (3.0, 3.0, 1.0))
        end

        @testset "Capsule" begin
            cap = FillableCapsule((0.0, 1.0, 0.0), (2.0, -1.0, 3.0), 0.5, 1.0)
            lo, hi = bounding_box(cap)
            @test lo == (-0.5, -1.5, -0.5)
            @test hi == (2.5, 1.5, 3.5)
            check_conservative(cap, (1.0, 0.0, 1.5), (2.0, 2.0, 3.0))
        end

        @testset "Slab (axis-aligned)" begin
            s = FillableSlab((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), 0.5, 1.0)
            lo, hi = bounding_box(s)
            @test lo == (-Inf, -Inf, 0.5)
            @test hi == (Inf, Inf, 1.5)
            check_conservative(s, (0.0, 0.0, 1.0), (5.0, 5.0, 1.0))
        end

        @testset "Slab (not axis-aligned)" begin
            s = FillableSlab((0.0, 0.0, 0.0), (1.0, 1.0, 0.0), 0.5, 1.0)
            lo, hi = bounding_box(s)
            @test all(isinf, lo)
            @test all(isinf, hi)
        end

        @testset "HalfSpace" begin
            # normal = +z: inside is z <= point z, bounded above only
            h = FillableHalfSpace((0.0, 0.0, 2.0), (0.0, 0.0, 1.0), 1.0)
            lo, hi = bounding_box(h)
            @test lo == (-Inf, -Inf, -Inf)
            @test hi == (Inf, Inf, 2.0)
            check_conservative(h, (0.0, 0.0, -3.0), (5.0, 5.0, 5.0))

            # normal = -z: inside is z >= point z, bounded below only
            h2 = FillableHalfSpace((0.0, 0.0, 2.0), (0.0, 0.0, -1.0), 1.0)
            lo2, hi2 = bounding_box(h2)
            @test lo2 == (-Inf, -Inf, 2.0)
            @test hi2 == (Inf, Inf, Inf)
            check_conservative(h2, (0.0, 0.0, 7.0), (5.0, 5.0, 5.0))

            # normal not axis-aligned: unbounded on every axis
            h3 = FillableHalfSpace((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 1.0)
            lo3, hi3 = bounding_box(h3)
            @test all(isinf, lo3)
            @test all(isinf, hi3)
        end
    end
end
