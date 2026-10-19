# GPU back-ends: rasterizing on a KernelAbstractions backend must produce the
# same voxel values as the CPU path, and the result must live on that backend.
# OpenCL runs on the pocl CPU device, so the GPU compile path is always tested.
# CUDA and Metal run when a device is available.

using pocl_jll, OpenCL
using CUDA
using Metal

backends = Backend[OpenCLBackend()]
CUDA.functional() && push!(backends, CUDABackend())
Metal.functional() && push!(backends, MetalBackend())

@testset "$(nameof(typeof(backend)))" for backend in backends
    # The shapes below are Float32, so every rasterize call needs
    # precision=Float32 to match.
    unitgrid() = Region((3, 3, 3), (1 // 1, 1 // 1, 1 // 1), (3 // 2, 3 // 2, 3 // 2))
    gpu(geometry, grid) = rasterize(geometry, grid, backend; precision=Float32)
    cpu(geometry, grid) = rasterize(geometry, grid; precision=Float32)

    @testset "rasterize on the backend matches the CPU for a sphere" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        arr = gpu(geometry, unitgrid())
        @test arr isa AbstractArray{Float32,3}
        @test get_backend(arr) == backend
        @test Array(arr) ≈ cpu(geometry, unitgrid())
    end

    @testset "rasterize! into a preallocated backend array matches the CPU" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        arr = KernelAbstractions.zeros(backend, Float32, 3, 3, 3)
        rasterize!(arr, geometry, unitgrid(); precision=Float32)
        @test Array(arr) ≈ cpu(geometry, unitgrid())
    end

    @testset "empty geometry is all background on the backend" begin
        geometry = Geometry((), 7.0f0)
        reg = Region((4, 4, 4), (1 // 1, 1 // 1, 1 // 1))
        @test all(==(7.0f0), Array(gpu(geometry, reg)))
    end

    @testset "shape ordering is preserved on the backend" begin
        c1 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 1.0f0)
        c2 = FillableCube((0.5f0, 0.5f0, 0.5f0), 10.0f0, 2.0f0)
        reg = Region((1, 1, 1), (1 // 1, 1 // 1, 1 // 1), (1 // 2, 1 // 2, 1 // 2))
        @test Array(gpu(Geometry([c1, c2], 0.0f0), reg))[1, 1, 1] == 1.0f0
        @test Array(gpu(Geometry([c2, c1], 0.0f0), reg))[1, 1, 1] == 2.0f0
    end

    @testset "anti-aliasing on the backend matches the CPU" begin
        sphere = FillableSphere((1.5f0, 1.5f0, 1.5f0), 1.0f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0, SuperResolutionAntiAliasing(2))
        host = Array(gpu(geometry, unitgrid()))
        @test host ≈ cpu(geometry, unitgrid())
        @test any(v -> 0.0f0 < v < 1.0f0, host)  # boundary blending occurred
    end

    @testset "multiple shapes and CSG on the backend match the CPU" begin
        a = FillableCube((1.5f0, 1.5f0, 1.5f0), 1.2f0, 1.0f0)
        b = FillableSphere((1.5f0, 1.5f0, 1.5f0), 0.6f0, 2.0f0)
        geometry = Geometry([csg_union(a, b)], 0.0f0)
        @test Array(gpu(geometry, unitgrid())) ≈ cpu(geometry, unitgrid())
    end

    @testset "composite grid on the backend matches the CPU" begin
        grid = refine(Region((8, 8, 8), (1 // 16, 1 // 16, 1 // 16)),
                      ((0 // 1, 0 // 1, 0 // 1), (1 // 8, 1 // 8, 1 // 8)))
        sphere = FillableSphere((0.0f0, 0.0f0, 0.0f0), 0.05f0, 1.0f0)
        geometry = Geometry([sphere], 0.0f0)
        field = gpu(geometry, grid)
        @test field isa CompositeField{Float32}
        @test get_backend(vec(field)) == backend
        @test length(field) == length(grid)
        host = cpu(geometry, grid)
        @test Array(vec(field)) ≈ vec(host)
        # regrid copies to the host first, so it works on a field on any backend.
        @test regrid(field) ≈ regrid(host)
    end
end
