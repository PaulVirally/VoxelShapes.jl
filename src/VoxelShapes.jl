module VoxelShapes

export add_shape

using KernelAbstractions
using CUDA

include("Types.jl")
using .Types
export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation
export interpolation, sdf, has_exact_sdf, center, bounding_box

include("Grids.jl")
using .Grids
export Region, CompositeGrid, refine, cells, scale, lower_corner, upper_corner, voxel_centers
export regions, nregions, finest, coordinates, cellvolumes

include("Interpolations.jl")
using .Interpolations
export LinearInterpolation, HarmonicInterpolation, GeometricMeanInterpolation, MaxInterpolation, MinInterpolation, DielectricInterpolation, MetalInterpolation
export interp_init, interp_accumulate, interp_finalize

include("Ellipsoids.jl")
using .Ellipsoids
export FillableEllipsoid, FillableSphere, radii

include("Cuboids.jl")
using .Cuboids
export FillableCuboid, FillableCube, half_lengths, lengths

include("Cylinders.jl")
using .Cylinders
export FillableCylinder, radius, half_height

include("Slabs.jl")
using .Slabs
export FillableSlab, FillableHalfSpace

include("Cones.jl")
using .Cones
export FillableCone

include("Tori.jl")
using .Tori
export FillableTorus, major_radius, minor_radius

include("Capsules.jl")
using .Capsules
export FillableCapsule

include("Transforms.jl")
using .Transforms
export Rotated

include("CSG.jl")
using .CSG
export UnionShape, IntersectionShape, DifferenceShape, ComplementShape
export csg_union, csg_intersect, csg_diff, csg_complement

include("Fills.jl")
using .Fills
export ConstantFill, RadialGradient, AxialGradient

include("AntiAliasing.jl")
using .AntiAliasing
export aa, NoAntiAliasing, SuperResolutionAntiAliasing, GaussianAntiAliasing
export SubpixelAntiAliasing, AdaptiveAntiAliasing

include("Geometry.jl")
export Geometry, sample

include("Rasterize.jl")
export rasterize, rasterize!, CompositeField, regionview, eachregion, regrid

end # module VoxelShapes
