# Interpolations

When anti-aliasing blends sub-voxel samples, the interpolation strategy controls how they are combined. Each shape carries its own, set at construction time via the `interpolation` keyword argument.

```julia
sphere = FillableSphere((0.5, 0.5, 0.5), 0.3, χ; interpolation=MetalInterpolation())
```

Interpolation only matters when anti-aliasing is active and a boundary voxel is being blended. For interior and exterior voxels, the fill value or background is returned directly.

## Strategies

`LinearInterpolation()` is the weighted arithmetic mean, and the default. Works for any scalar or vector fill value.

`HarmonicInterpolation()` is the weighted harmonic mean. Use it for positive-definite quantities. A single zero sub-sample collapses the voxel to zero.

`GeometricMeanInterpolation()` is the weighted geometric mean, computed as `exp(∑ wᵢ log(vᵢ))`. Positive-definite quantities only.

`MaxInterpolation()` takes the largest sub-sample and ignores weights. If any sub-sample lands inside the shape, the voxel takes the shape's fill value. Good when partial presence should count as full presence.

`MinInterpolation()` takes the smallest sub-sample. The voxel only gets the fill value when every sub-sample is inside. The conservative choice: the shape has to fully enclose the voxel before it claims it.

`DielectricInterpolation()` is linear interpolation of electric susceptibility χ. Identical to `LinearInterpolation`, just named for clarity in electromagnetic work.

`MetalInterpolation()` handles metallic media. It interpolates the complex refractive index ñ = sqrt(1 + χ), then recovers χ = ñ² - 1. From Christiansen et al. (2019). Pass χ as the fill value, not ñ.

## Summary table

| Strategy | Aggregation | When to use |
|---|---|---|
| `LinearInterpolation` | Arithmetic mean | General purpose |
| `HarmonicInterpolation` | Harmonic mean | Positive quantities with strong variation |
| `GeometricMeanInterpolation` | Geometric mean | Positive quantities spanning orders of magnitude |
| `MaxInterpolation` | Max | Inclusive: any sub-sample inside claims the voxel |
| `MinInterpolation` | Min | Conservative: all sub-samples must be inside |
| `DielectricInterpolation` | Arithmetic mean | Electric susceptibility (Maxwell dielectrics) |
| `MetalInterpolation` | Via ñ = sqrt(1+χ) | Complex susceptibility of metals |

See the [API reference](@ref "API reference") for full signatures.
