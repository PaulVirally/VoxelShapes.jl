# CUDA back-end: CuArray(world) and fill!(::CuArray, world) must produce the
# same voxel values as the CPU path, and the returned array must live on the GPU.
# All tests are skipped when no CUDA device is available.

using CUDA

if !CUDA.functional()
    @warn "No CUDA device found: skipping CUDA tests"
else

@testset "CUDA back-end" begin

    @testset "CuArray(world) matches Array(world) for a sphere" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0f0)
        w = World((3, 3, 3), (1.0f0, 1.0f0, 1.0f0), [sphere], 0.0f0, NoAntiAliasing())
        cpu = Array(w)
        gpu = CuArray(w)
        @test gpu isa CuArray{Float32, 3}
        @test Array(gpu) ≈ cpu
    end

    @testset "fill! into a preallocated CuArray matches Array(world)" begin
        sphere = FillableSphere((1.5, 1.5, 1.5), 0.6, 1.0f0)
        w = World((3, 3, 3), (1.0f0, 1.0f0, 1.0f0), [sphere], 0.0f0, NoAntiAliasing())
        arr = CUDA.zeros(Float32, 3, 3, 3)
        fill!(arr, w)
        @test Array(arr) ≈ Array(w)
    end

    @testset "empty world is all background on GPU" begin
        w = World((4, 4, 4), (1.0f0, 1.0f0, 1.0f0), (), 7.0f0, NoAntiAliasing())
        @test all(==(7.0f0), Array(CuArray(w)))
    end

    @testset "shape ordering is preserved on GPU" begin
        c1 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 1.0f0)
        c2 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 2.0f0)
        w12 = World((1, 1, 1), (1.0f0, 1.0f0, 1.0f0), [c1, c2], 0.0f0, NoAntiAliasing())
        w21 = World((1, 1, 1), (1.0f0, 1.0f0, 1.0f0), [c2, c1], 0.0f0, NoAntiAliasing())
        @test Array(CuArray(w12))[1, 1, 1] == 1.0f0
        @test Array(CuArray(w21))[1, 1, 1] == 2.0f0
    end

    @testset "CuArray(world) with SuperResolutionAntiAliasing matches CPU" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 1.0f0, 1.0f0)
        w = World((3, 3, 3), (1.0f0, 1.0f0, 1.0f0), [sphere], 0.0f0, SuperResolutionAntiAliasing(2))
        cpu = Array(w)
        gpu = Array(CuArray(w))
        @test gpu ≈ cpu
        @test any(v -> 0.0f0 < v < 1.0f0, gpu)  # boundary blending occurred
    end

    @testset "CuArray(world) with multiple shapes and CSG matches CPU" begin
        a = FillableCube((1.5f0, 1.5f0, 1.5f0), 1.2f0, 1.0f0)
        b = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 2.0f0)
        u = csg_union(a, b)
        w = World((3, 3, 3), (1.0f0, 1.0f0, 1.0f0), [u], 0.0f0, NoAntiAliasing())
        @test Array(CuArray(w)) ≈ Array(w)
    end

end # @testset "CUDA back-end"

end # if CUDA.functional()
