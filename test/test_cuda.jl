# CUDA back-end: rasterizing onto a CuArray must produce the same voxel values
# as the CPU path, and the result must live on the GPU. All tests are skipped
# when no CUDA device is available.

using CUDA

if !CUDA.functional()
    @warn "No CUDA device found: skipping CUDA tests"
else

@testset "CUDA back-end" begin
    # The shapes below are Float32, so every rasterize call needs
    # precision=Float32 to match.
    unitgrid() = Region((3, 3, 3), (1 // 1, 1 // 1, 1 // 1), (3 // 2, 3 // 2, 3 // 2))
    gpu(geometry, grid) = rasterize(geometry, grid, CuArray; precision=Float32)
    cpu(geometry, grid) = rasterize(geometry, grid; precision=Float32)

    @testset "rasterize onto a CuArray matches the CPU for a sphere" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        arr = gpu(geometry, unitgrid())
        @test arr isa CuArray{Float32,3}
        @test Array(arr) ≈ cpu(geometry, unitgrid())
    end

    @testset "rasterize! into a preallocated CuArray matches the CPU" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        arr = CUDA.zeros(Float32, 3, 3, 3)
        rasterize!(arr, geometry, unitgrid(); precision=Float32)
        @test Array(arr) ≈ cpu(geometry, unitgrid())
    end

    @testset "empty geometry is all background on the GPU" begin
        geometry = Geometry((), 7.0f0)
        reg = Region((4, 4, 4), (1 // 1, 1 // 1, 1 // 1))
        @test all(==(7.0f0), Array(gpu(geometry, reg)))
    end

    @testset "shape ordering is preserved on the GPU" begin
        c1 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 1.0f0)
        c2 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 2.0f0)
        reg = Region((1, 1, 1), (1 // 1, 1 // 1, 1 // 1), (1 // 2, 1 // 2, 1 // 2))
        @test Array(gpu(Geometry([c1, c2], 0.0f0), reg))[1, 1, 1] == 1.0f0
        @test Array(gpu(Geometry([c2, c1], 0.0f0), reg))[1, 1, 1] == 2.0f0
    end

    @testset "anti-aliasing on the GPU matches the CPU" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 1.0f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0, SuperResolutionAntiAliasing(2))
        host = Array(gpu(geometry, unitgrid()))
        @test host ≈ cpu(geometry, unitgrid())
        @test any(v -> 0.0f0 < v < 1.0f0, host)  # boundary blending occurred
    end

    @testset "multiple shapes and CSG on the GPU match the CPU" begin
        a = FillableCube((1.5f0, 1.5f0, 1.5f0), 1.2f0, 1.0f0)
        b = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 2.0f0)
        geometry = Geometry([csg_union(a, b)], 0.0f0)
        @test Array(gpu(geometry, unitgrid())) ≈ cpu(geometry, unitgrid())
    end

    @testset "composite grid on the GPU matches the CPU" begin
        grid = refine(Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16)),
                      ((0 // 1, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8)))
        sphere = FillableSphere((0.0f0, 0.0f0, 0.0f0), 0.05f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        field = gpu(geometry, grid)
        @test field isa CompositeField{Float32}
        @test vec(field) isa CuVector{Float32}
        @test length(field) == length(grid)
        host = cpu(geometry, grid)
        @test Array(vec(field)) ≈ vec(host)
        # regrid copies to the host first, so it works on a GPU field.
        @test regrid(field) ≈ regrid(host)
    end
end # @testset "CUDA back-end"

end # if CUDA.functional()
