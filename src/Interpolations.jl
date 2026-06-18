module Interpolations

using ..Types

export LinearInterpolation
export interp_init, interp_accumulate, interp_finalize

# Fold contract
#
# interp_init(interp, ::Type{U}) -> initial accumulator of type U
# interp_accumulate(interp, acc, val, w) -> updated accumulator
# interp_finalize(interp, acc) -> final mixed value
#
# The sampler calls this trio once per voxel, folding (value, weight) pairs
# where weights sum to 1 and value is the per-sub-sample fill-or-background.

# LinearInterpolation
struct LinearInterpolation <: AbstractInterpolation end

interp_init(::LinearInterpolation, ::Type{U}) where {U} = zero(U)
interp_accumulate(::LinearInterpolation, acc, val, weight) = acc + weight * val
interp_finalize(::LinearInterpolation, acc) = acc

end # module Interpolations
