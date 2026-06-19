module Fills

export ConstantFill, RadialGradient, AxialGradient

"""
    ConstantFill{V}

Fill function that returns the same value everywhere inside a shape.

Named `isbits`-compatible alternative to the `_ -> fill_val` closure, which
avoids heap allocation on GPUs.

# Fields
- `value`: the constant fill value
"""
struct ConstantFill{V}
    value::V
end
(f::ConstantFill)(_) = f.value

"""
    RadialGradient{V}

Fill function that interpolates linearly from `inner_value` at the shape center
to `outer_value` at the shape surface.

The radial coordinate is `r = norm(local_coords)`, clamped to `[0, 1]`.
Local coordinates follow each shape's own convention: for spheres and ellipsoids,
`r = 1` at the surface; for cylinders the radial and axial coordinates are
independent, so behavior at `r > 1` may be unexpected.

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

Fill function that interpolates linearly along one local axis of a shape.

The axis coordinate ranges from -1 to +1 across the shape. `v0` is returned
at `-1` and `v1` at `+1`; values are clamped outside that range.

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
