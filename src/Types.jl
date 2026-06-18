module Types

export AbstractFillableShape, AbstractAntiAliasing, AbstractInterpolation, interpolation

abstract type AbstractFillableShape end
abstract type AbstractAntiAliasing end

"""
    abstract type AbstractInterpolation

Abstract supertype for all per-sample fold (init / accumulate / finalize) mixing rules.
Implement `interp_init`, `interp_accumulate`, and `interp_finalize` for each concrete subtype.
"""
abstract type AbstractInterpolation end

"""
    interpolation(shape) -> AbstractInterpolation

Return the interpolation scheme attached to `shape`.
"""
function interpolation end

end
