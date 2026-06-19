using Test
using Aqua
using VoxelShapes
using StaticArrays
using LinearAlgebra: norm

# The suite is split into two passes.
#
# Pass 1 (TDD): tests written against the documented behavior of the public
# API, treating the implementation as a black box. Grouped by feature area.
#
# Pass 2 (integration): tests added after reading the implementation, covering
# interactions between components that the black-box tests did not reach
# (anti-aliasing x interpolation, CSG x SDF combinators, rotation isometry,
# the World rasterization pipeline, etc.).

@testset "VoxelShapes.jl" begin
    include("test_aqua.jl")

    # --- Pass 1: feature-by-feature ---
    include("test_interpolations.jl")
    include("test_shapes.jl")
    include("test_fills.jl")
    include("test_transforms.jl")
    include("test_csg.jl")
    include("test_antialiasing.jl")
    include("test_world.jl")

    # --- Pass 2: cross-cutting integration ---
    include("test_integration.jl")

    # --- GPU ---
    include("test_cuda.jl")
end
