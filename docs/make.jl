using Documenter
using VoxelShapes

makedocs(
    sitename = "VoxelShapes.jl",
    repo = "github.com/PaulVirally/VoxelShapes.jl.git",
    modules  = [
        VoxelShapes,
        VoxelShapes.Types,
        VoxelShapes.Interpolations,
        VoxelShapes.Ellipsoids,
        VoxelShapes.Cuboids,
        VoxelShapes.Cylinders,
        VoxelShapes.Slabs,
        VoxelShapes.Cones,
        VoxelShapes.Tori,
        VoxelShapes.Capsules,
        VoxelShapes.Transforms,
        VoxelShapes.CSG,
        VoxelShapes.Fills,
        VoxelShapes.AntiAliasing,
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://paulvirally.github.io/VoxelShapes.jl",
    ),
    remotes  = nothing,
    pages = [
        "Getting started"  => "index.md",
        "Shapes"           => "shapes.md",
        "Fill functions"   => "fills.md",
        "Anti-aliasing"    => "antialiasing.md",
        "Interpolations"   => "interpolations.md",
        "CSG"              => "csg.md",
        "Rotation"         => "transforms.md",
        "GPU"              => "gpu.md",
        "Extending"        => "extending.md",
        "API reference"    => "api.md",
    ],
    checkdocs = :none,
    warnonly  = true,
)

deploydocs(
    repo       = "github.com/PaulVirally/VoxelShapes.jl",
    devbranch  = "main",
)
