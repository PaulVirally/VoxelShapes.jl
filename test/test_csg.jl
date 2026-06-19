# Constructive solid geometry: union, intersection, difference, complement.
# Containment follows boolean logic; SDFs use the classic min/max combinators;
# fill delegates to the primary operand.

@testset "CSG" begin
    # Two overlapping unit spheres on the x-axis, distinct fill values.
    A = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)
    B = FillableSphere((1.5, 0.0, 0.0), 1.0, 2.0)

    @testset "constructors return the right wrapper types" begin
        @test csg_union(A, B) isa UnionShape
        @test csg_intersect(A, B) isa IntersectionShape
        @test csg_diff(A, B) isa DifferenceShape
        @test csg_complement(A) isa ComplementShape
        @test csg_union(A, B) isa AbstractFillableShape
    end

    @testset "union: inside a or b" begin
        u = csg_union(A, B)
        @test (0.0, 0.0, 0.0) in u        # in A
        @test (1.5, 0.0, 0.0) in u        # in B
        @test !((5.0, 0.0, 0.0) in u)     # in neither
        # fill delegates to A where the point is in A, else to B
        @test fill(u, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)) == 1.0
        @test fill(u, (1.5, 0.0, 0.0), (1.0, 1.0, 1.0)) == 2.0   # only in B
    end

    @testset "intersection: inside both" begin
        i = csg_intersect(A, B)
        @test (0.75, 0.0, 0.0) in i        # in both (overlap region)
        @test !((0.0, 0.0, 0.0) in i)      # in A only
        @test !((1.5, 0.0, 0.0) in i)      # in B only
        @test fill(i, (0.75, 0.0, 0.0), (1.0, 1.0, 1.0)) == 1.0  # delegates to A
    end

    @testset "difference: in a but not b" begin
        d = csg_diff(A, B)
        @test (0.0, 0.0, 0.0) in d         # in A, not in B
        @test !((0.75, 0.0, 0.0) in d)     # in both -> carved away
        @test !((1.5, 0.0, 0.0) in d)      # in B only
        @test fill(d, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)) == 1.0
    end

    @testset "complement: everything outside" begin
        c = csg_complement(A)
        @test (5.0, 0.0, 0.0) in c
        @test !((0.0, 0.0, 0.0) in c)
    end

    @testset "SDF combinators (using exact-SDF cubes)" begin
        # Two unit cubes (half-length 1) for exact SDFs.
        P = FillableCube((0.0, 0.0, 0.0), 2.0, 1.0)
        Q = FillableCube((3.0, 0.0, 0.0), 2.0, 1.0)
        @test sdf(csg_union(P, Q), (5.0, 0.0, 0.0)) ≈ min(sdf(P, (5.0,0.0,0.0)), sdf(Q, (5.0,0.0,0.0)))
        @test sdf(csg_intersect(P, Q), (1.5, 0.0, 0.0)) ≈ max(sdf(P, (1.5,0.0,0.0)), sdf(Q, (1.5,0.0,0.0)))
        @test sdf(csg_diff(P, Q), (0.0, 0.0, 0.0)) ≈ max(sdf(P, (0.0,0.0,0.0)), -sdf(Q, (0.0,0.0,0.0)))
        @test sdf(csg_complement(P), (5.0, 0.0, 0.0)) ≈ -sdf(P, (5.0, 0.0, 0.0))
    end

    @testset "has_exact_sdf requires all operands exact" begin
        P = FillableCube((0.0, 0.0, 0.0), 2.0, 1.0)         # exact
        Q = FillableCube((3.0, 0.0, 0.0), 2.0, 1.0)         # exact
        sph = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0)     # not exact
        @test has_exact_sdf(csg_union(P, Q)) == true
        @test has_exact_sdf(csg_intersect(P, Q)) == true
        @test has_exact_sdf(csg_diff(P, Q)) == true
        @test has_exact_sdf(csg_complement(P)) == true
        @test has_exact_sdf(csg_union(P, sph)) == false
        @test has_exact_sdf(csg_complement(sph)) == false
    end

    @testset "interpolation delegates to the primary operand" begin
        A2 = FillableSphere((0.0, 0.0, 0.0), 1.0, 1.0; interpolation=MaxInterpolation())
        @test interpolation(csg_union(A2, B)) isa MaxInterpolation
        @test interpolation(csg_intersect(A2, B)) isa MaxInterpolation
        @test interpolation(csg_diff(A2, B)) isa MaxInterpolation
        @test interpolation(csg_complement(A2)) isa MaxInterpolation
    end

    @testset "nesting" begin
        # (A ∪ B) \ C should still behave sensibly.
        C = FillableSphere((0.0, 0.0, 0.0), 0.5, 3.0)
        nested = csg_diff(csg_union(A, B), C)
        @test !((0.0, 0.0, 0.0) in nested)   # carved out by C
        @test (1.5, 0.0, 0.0) in nested      # in B, outside C
    end
end
