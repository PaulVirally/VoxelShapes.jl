# Interpolation fold protocol: interp_init / interp_accumulate / interp_finalize.
# Each interpolation combines weighted (value, weight) sub-voxel samples into one
# blended value. Weights across all samples sum to one.

# Helper: fold a list of (value, weight) pairs through an interpolation scheme.
function blend(interp, ::Type{U}, samples) where {U}
    acc = interp_init(interp, U)
    for (val, w) in samples
        acc = interp_accumulate(interp, acc, val, w)
    end
    return interp_finalize(interp, acc)
end

@testset "Interpolations" begin
    @testset "types" begin
        for I in (LinearInterpolation, HarmonicInterpolation, GeometricMeanInterpolation,
                  MaxInterpolation, MinInterpolation, DielectricInterpolation, MetalInterpolation)
            @test I() isa AbstractInterpolation
        end
    end

    @testset "LinearInterpolation = weighted arithmetic mean" begin
        li = LinearInterpolation()
        @test blend(li, Float64, ((1.0, 0.5), (3.0, 0.5))) ≈ 2.0
        @test blend(li, Float64, ((10.0, 0.25), (2.0, 0.75))) ≈ 4.0
        # A single full-weight sample returns that sample unchanged.
        @test blend(li, Float64, ((7.5, 1.0),)) ≈ 7.5
        # The initial accumulator is the additive identity.
        @test interp_init(li, Float64) == 0.0
    end

    @testset "HarmonicInterpolation = weighted harmonic mean" begin
        hi = HarmonicInterpolation()
        # harmonic mean of 1 and 3 with equal weight = 1.5
        @test blend(hi, Float64, ((1.0, 0.5), (3.0, 0.5))) ≈ 1.5
        @test blend(hi, Float64, ((4.0, 1.0),)) ≈ 4.0
        # A zero sub-sample collapses the result to zero (documented behavior).
        @test blend(hi, Float64, ((0.0, 0.5), (2.0, 0.5))) == 0.0
    end

    @testset "GeometricMeanInterpolation = weighted geometric mean" begin
        gi = GeometricMeanInterpolation()
        # geometric mean of 1 and 4 with equal weight = 2
        @test blend(gi, Float64, ((1.0, 0.5), (4.0, 0.5))) ≈ 2.0
        @test blend(gi, Float64, ((8.0, 1.0),)) ≈ 8.0
        # geometric mean of 2 and 8 = 4
        @test blend(gi, Float64, ((2.0, 0.5), (8.0, 0.5))) ≈ 4.0
    end

    @testset "MaxInterpolation ignores weights" begin
        mi = MaxInterpolation()
        @test blend(mi, Float64, ((1.0, 0.9), (3.0, 0.1))) ≈ 3.0
        @test blend(mi, Float64, ((5.0, 0.5), (2.0, 0.5))) ≈ 5.0
    end

    @testset "MinInterpolation ignores weights" begin
        mi = MinInterpolation()
        @test blend(mi, Float64, ((1.0, 0.1), (3.0, 0.9))) ≈ 1.0
        @test blend(mi, Float64, ((5.0, 0.5), (2.0, 0.5))) ≈ 2.0
    end

    @testset "DielectricInterpolation behaves like Linear" begin
        di = DielectricInterpolation()
        @test blend(di, Float64, ((1.0, 0.5), (3.0, 0.5))) ≈ 2.0
        @test blend(di, Float64, ((10.0, 0.25), (2.0, 0.75))) ≈ 4.0
    end

    @testset "MetalInterpolation via complex refractive index" begin
        # Fill value is χ; ñ = sqrt(1 + χ); interpolate ñ; recover χ = ñ² - 1.
        mp = MetalInterpolation()
        # A single full-weight sample is the identity: χ in, χ out.
        @test blend(mp, Float64, ((3.0, 1.0),)) ≈ 3.0
        # Two identical χ=3 samples (ñ=2) also recover χ=3.
        @test blend(mp, Float64, ((3.0, 0.5), (3.0, 0.5))) ≈ 3.0
        # Distinct samples: ñ interpolation of χ=0 (ñ=1) and χ=3 (ñ=2),
        # equal weight -> ñ=1.5 -> χ = 1.5² - 1 = 1.25
        @test real(blend(mp, Float64, ((0.0, 0.5), (3.0, 0.5)))) ≈ 1.25
    end
end
