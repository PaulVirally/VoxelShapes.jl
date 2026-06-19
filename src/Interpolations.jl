module Interpolations

using ..Types

export LinearInterpolation
export HarmonicInterpolation, GeometricMeanInterpolation, MaxInterpolation, MinInterpolation, DielectricInterpolation, MetalInterpolation
export interp_init, interp_accumulate, interp_finalize

"""
    interp_init(interp, ::Type{U}) -> acc

Return the initial accumulator value for an interpolation fold over sub-voxel samples.

`U` is the element type of the fill value.
"""
function interp_init end

"""
    interp_accumulate(interp, acc, val, weight) -> acc

Fold one `(val, weight)` pair into `acc`.

Called once per sub-voxel sample. Weights across all samples sum to one.
`val` is either the shape's fill value or the background, depending on containment.
"""
function interp_accumulate end

"""
    interp_finalize(interp, acc) -> result

Convert the accumulated state into the final blended voxel value.
"""
function interp_finalize end

"""
    LinearInterpolation <: AbstractInterpolation

Weighted arithmetic mean of sub-voxel samples. The standard choice for scalar
and vector fill values.
"""
struct LinearInterpolation <: AbstractInterpolation end

interp_init(::LinearInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::LinearInterpolation, acc, val, weight) = acc + weight * val
interp_finalize(::LinearInterpolation, acc) = acc

"""
    HarmonicInterpolation <: AbstractInterpolation

Weighted harmonic mean of sub-voxel samples.

Intended for strictly positive quantities. A zero sub-sample value collapses
the voxel to zero (by design, via the `w/val → Inf → 1/Inf = 0` path).
"""
struct HarmonicInterpolation <: AbstractInterpolation end

interp_init(::HarmonicInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::HarmonicInterpolation, acc, val, w) = acc + w / val # val==0 -> Inf
interp_finalize(::HarmonicInterpolation, acc) = one(typeof(acc)) / acc # Inf -> 0

"""
    GeometricMeanInterpolation <: AbstractInterpolation

Weighted geometric mean of sub-voxel samples.

Intended for strictly positive quantities only. Implemented via
`exp(∑ wᵢ log(vᵢ))`.
"""
struct GeometricMeanInterpolation <: AbstractInterpolation end

interp_init(::GeometricMeanInterpolation, ::Type{U}) where {U} = zero(float(U))
interp_accumulate(::GeometricMeanInterpolation, acc, val, w) = acc + w * log(val)
interp_finalize(::GeometricMeanInterpolation, acc) = exp(acc)

"""
    MaxInterpolation <: AbstractInterpolation

Take the maximum sub-voxel sample, ignoring weights.

Useful when any sub-sample belonging to a shape should claim the whole voxel
(for example: if any sub-sample is metal, the voxel should be metal).
"""
struct MaxInterpolation <: AbstractInterpolation end

interp_init(::MaxInterpolation, ::Type{U}) where {U} = typemin(U)
interp_accumulate(::MaxInterpolation, acc, val, _) = max(acc, val)
interp_finalize(::MaxInterpolation, acc) = acc

"""
    MinInterpolation <: AbstractInterpolation

Take the minimum sub-voxel sample, ignoring weights.

Conservative inclusion: a voxel is assigned a fill value only when all
sub-samples fall inside the shape.
"""
struct MinInterpolation <: AbstractInterpolation end

interp_init(::MinInterpolation, ::Type{U}) where {U} = typemax(U)
interp_accumulate(::MinInterpolation, acc, val, _) = min(acc, val)
interp_finalize(::MinInterpolation, acc) = acc

"""
    DielectricInterpolation <: AbstractInterpolation

Linear interpolation of electric susceptibility χ for Maxwell dielectric filling.

Equivalent to `LinearInterpolation` but named for clarity in electromagnetic
simulation workflows where the fill value represents χ.
"""
struct DielectricInterpolation <: AbstractInterpolation end
interp_init(::DielectricInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::DielectricInterpolation, acc, val, w) = acc + w * val
interp_finalize(::DielectricInterpolation, acc) = acc

"""
    MetalInterpolation <: AbstractInterpolation

Nonlinear interpolation of complex electric susceptibility χ for metallic media.

Linearly interpolates the complex refractive index ñ = n + iκ, then recovers
χ = ñ² - 1. This gives Re(χ) = n² - κ² - 1 and Im(χ) = 2nκ.

The fill value should be χ (not ñ). Method from Christiansen et al. (2019).
"""
struct MetalInterpolation <: AbstractInterpolation end
interp_init(::MetalInterpolation, ::Type{U}) where {U} = zero(complex(float(U)))
interp_accumulate(::MetalInterpolation, acc, val, w) = acc + w * sqrt(complex(one(val) + val))
interp_finalize(::MetalInterpolation, acc) = acc^2 - 1

end # module Interpolations
