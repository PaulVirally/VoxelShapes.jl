module Interpolations

using ..Types

export LinearInterpolation
export HarmonicInterpolation, GeometricMeanInterpolation, MaxInterpolation, MinInterpolation, DielectricInterpolation, MetalInterpolation
export interp_init, interp_accumulate, interp_finalize

"""
    interp_init(interp, ::Type{U}) -> acc

Initial accumulator value for a fold over sub-voxel samples.

# Arguments
- `U`: element type of the fill value
"""
function interp_init end

"""
    interp_accumulate(interp, acc, val, weight) -> acc

Fold one `(val, weight)` pair into `acc`, once per sub-voxel sample.

# Arguments
- `val`: the shape's fill value or the background, depending on containment
- `weight`: weights across all samples sum to one
"""
function interp_accumulate end

"""
    interp_finalize(interp, acc) -> result

Convert the accumulated state into the final blended voxel value.
"""
function interp_finalize end

"""
    LinearInterpolation <: AbstractInterpolation

Weighted arithmetic mean of sub-voxel samples.
"""
struct LinearInterpolation <: AbstractInterpolation end

interp_init(::LinearInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::LinearInterpolation, acc, val, weight) = acc + weight * val
interp_finalize(::LinearInterpolation, acc) = acc

"""
    HarmonicInterpolation <: AbstractInterpolation

Weighted harmonic mean of sub-voxel samples, for strictly positive quantities.

A zero sub-sample collapses the voxel to zero, via `w/val → Inf → 1/Inf = 0`.
"""
struct HarmonicInterpolation <: AbstractInterpolation end

interp_init(::HarmonicInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::HarmonicInterpolation, acc, val, w) = acc + w / val # val==0 -> Inf
interp_finalize(::HarmonicInterpolation, acc) = one(typeof(acc)) / acc # Inf -> 0

"""
    GeometricMeanInterpolation <: AbstractInterpolation

Weighted geometric mean of sub-voxel samples, for strictly positive
quantities, via `exp(∑ wᵢ log(vᵢ))`.
"""
struct GeometricMeanInterpolation <: AbstractInterpolation end

interp_init(::GeometricMeanInterpolation, ::Type{U}) where {U} = zero(float(U))
interp_accumulate(::GeometricMeanInterpolation, acc, val, w) = acc + w * log(val)
interp_finalize(::GeometricMeanInterpolation, acc) = exp(acc)

"""
    MaxInterpolation <: AbstractInterpolation

Maximum sub-voxel sample, ignoring weights.
"""
struct MaxInterpolation <: AbstractInterpolation end

interp_init(::MaxInterpolation, ::Type{U}) where {U} = typemin(U)
interp_accumulate(::MaxInterpolation, acc, val, _) = max(acc, val)
interp_finalize(::MaxInterpolation, acc) = acc

"""
    MinInterpolation <: AbstractInterpolation

Minimum sub-voxel sample, ignoring weights.

A voxel gets a fill value only when all sub-samples are inside the shape.
"""
struct MinInterpolation <: AbstractInterpolation end

interp_init(::MinInterpolation, ::Type{U}) where {U} = typemax(U)
interp_accumulate(::MinInterpolation, acc, val, _) = min(acc, val)
interp_finalize(::MinInterpolation, acc) = acc

"""
    DielectricInterpolation <: AbstractInterpolation

Linear interpolation of electric susceptibility χ, equivalent to
`LinearInterpolation`.
"""
struct DielectricInterpolation <: AbstractInterpolation end
interp_init(::DielectricInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::DielectricInterpolation, acc, val, w) = acc + w * val
interp_finalize(::DielectricInterpolation, acc) = acc

"""
    MetalInterpolation <: AbstractInterpolation

Nonlinear interpolation of complex electric susceptibility χ for metals.

Linearly interpolates the complex refractive index ñ = n + iκ, then recovers
χ = ñ² - 1, giving Re(χ) = n² - κ² - 1 and Im(χ) = 2nκ. The fill value must
be χ, not ñ. Method from Christiansen et al. (2019).
"""
struct MetalInterpolation <: AbstractInterpolation end
interp_init(::MetalInterpolation, ::Type{U}) where {U} = zero(complex(float(U)))
interp_accumulate(::MetalInterpolation, acc, val, w) = acc + w * sqrt(complex(one(val) + val))
interp_finalize(::MetalInterpolation, acc) = acc^2 - 1

end # module Interpolations
