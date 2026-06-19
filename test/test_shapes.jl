# Per-primitive behavior: containment (`in`), fill value (`fill`), accessors,
# signed distance (`sdf`) and the `has_exact_sdf` flag. Geometry conventions are
# taken from the docstrings / README, not the implementation.

const VS = 1.0  # nominal voxel size used where fill needs a voxel_size argument

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
        c = FillableCube((0.0, 0.0, 0.0), 2.0, 9.0)   # side length 2 -> half-length 1
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
        @test sdf(c, (0.0, 0.0, 0.0)) ≈ -1.0     # 1 unit inside (to nearest face)
        @test sdf(c, (1.0, 0.0, 0.0)) ≈ 0.0      # on the surface
    end

    @testset "FillableCylinder" begin
        cyl = FillableCylinder((0.0, 0.0, 0.0), 1.0, 2.0, 3.0)  # axis=3 (z), fill 3
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
        @test sdf(cyl, (2.0, 0.0, 0.0)) ≈ 1.0    # 1 unit outside curved surface
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
        s = FillableSlab((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 1.0, 7.0)  # half_thickness 1
        @test (0.0, 0.0, 0.0) in s
        @test (0.0, 0.0, 1.0) in s               # on the slab face
        @test !((0.0, 0.0, 1.5) in s)
        @test (0.0, 0.0, -0.5) in s
        @test fill(s, (0.0, 0.0, 0.0), (VS, VS, VS)) == 7.0

        @test has_exact_sdf(s) == true
        @test sdf(s, (0.0, 0.0, 2.0)) ≈ 1.0      # 1 unit outside the slab
        @test sdf(s, (0.0, 0.0, 0.0)) ≈ -1.0     # center, 1 unit from each face

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
        @test sdf(cap, (1.0, 0.0, 1.0)) ≈ 0.5    # distance 1 to axis, minus radius
        @test sdf(cap, (0.0, 0.0, 1.0)) ≈ -0.5
    end
end
