# API reference

## World

```@docs
World
add_shape
```

## Shapes

```@docs
FillableSphere
FillableEllipsoid
radii
FillableCuboid
FillableCube
half_lengths
lengths
FillableCylinder
radius
half_height
FillableTorus
major_radius
minor_radius
FillableCapsule
FillableCone
FillableSlab
FillableHalfSpace
```

## Fill functions

```@docs
ConstantFill
RadialGradient
AxialGradient
```

## Anti-aliasing

```@docs
NoAntiAliasing
SuperResolutionAntiAliasing
SubpixelAntiAliasing
GaussianAntiAliasing
AdaptiveAntiAliasing
aa
```

## Interpolations

```@docs
LinearInterpolation
HarmonicInterpolation
GeometricMeanInterpolation
MaxInterpolation
MinInterpolation
DielectricInterpolation
MetalInterpolation
interp_init
interp_accumulate
interp_finalize
```

## CSG

```@docs
UnionShape
IntersectionShape
DifferenceShape
ComplementShape
csg_union
csg_intersect
csg_diff
csg_complement
```

## Transforms

```@docs
Rotated
```

## Abstract types and interface

```@docs
AbstractFillableShape
AbstractAntiAliasing
AbstractInterpolation
interpolation
sdf
has_exact_sdf
center
```
