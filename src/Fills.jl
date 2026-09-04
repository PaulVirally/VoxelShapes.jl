module Fills

export ConstantFill, RadialGradient, AxialGradient

"""
    ConstantFill{V}

`isbits`-compatible fill function that returns `value` everywhere.

# Fields
- `value`: the constant fill value
"""
struct ConstantFill{V}
    value::V
end
(f::ConstantFill)(_) = f.value

"""
    RadialGradient{V}

Fill function that interpolates linearly along `r = norm(local_coords)`,
clamped to `[0, 1]`.

# Fields
- `inner_value`: value at `r = 0`
- `outer_value`: value at `r = 1`
"""
struct RadialGradient{V}
    inner_value::V
    outer_value::V
end
function (f::RadialGradient)(local_coords::NTuple{3})
    r = clamp(sqrt(local_coords[1]^2 + local_coords[2]^2 + local_coords[3]^2), 0, 1)
    return f.inner_value + (f.outer_value - f.inner_value) * r
end

"""
    AxialGradient{V}

Fill function that interpolates linearly along one local axis of a shape,
clamped outside `[-1, 1]`.

# Fields
- `axis`: which local coordinate to interpolate along (1, 2, or 3)
- `v0`: value at local coordinate -1
- `v1`: value at local coordinate +1
"""
struct AxialGradient{V}
    axis::Int
    v0::V
    v1::V
end
function (f::AxialGradient)(local_coords::NTuple{3})
    t = clamp((local_coords[f.axis] + 1) / 2, 0, 1)
    return f.v0 + (f.v1 - f.v0) * t
end

end # module
