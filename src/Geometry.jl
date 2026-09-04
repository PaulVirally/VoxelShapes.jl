"""
    Geometry{U, S, A}

Immutable description of what to draw: a tuple of shapes, a background value,
and an anti-aliasing strategy.

Hand it to [`rasterize`](@ref) with a [`Region`](@ref) or
[`CompositeGrid`](@ref) to get an array, or call [`sample`](@ref) at a single
position. Shapes are evaluated in order. The first one whose containment test
succeeds claims the position. A shape that evaluates to `background` is
transparent: it falls through to the next shape. [`add_shape`](@ref) appends
a shape and returns a new `Geometry`.

# Fields
- `shapes::S`: tuple of `AbstractFillableShape` values, evaluated first to last
- `background::U`: value where no shape claims the position
- `aa::A`: anti-aliasing strategy applied at every sample
"""
struct Geometry{U, S<:Tuple, A<:AbstractAntiAliasing}
    shapes::S
    background::U
    aa::A
end

"""
    Geometry(shapes, background, aa=NoAntiAliasing())

Construct a geometry from a `Tuple` or `AbstractVector` of shapes.
"""
Geometry(shapes::Tuple, background) = Geometry(shapes, background, NoAntiAliasing())
Geometry(shapes::AbstractVector, background, antialiasing::AbstractAntiAliasing=NoAntiAliasing()) =
    Geometry(Tuple(shapes), background, antialiasing)

Base.eltype(geometry::Geometry) = typeof(geometry.background)

"""
    add_shape(geometry, shape) -> Geometry

Return a new `Geometry` with `shape` appended to the end of the shape tuple.
"""
add_shape(geometry::Geometry, shape) = Geometry((geometry.shapes..., shape), geometry.background, geometry.aa)

"""
    sample(geometry, pos::NTuple{3,T}, voxel_size::NTuple{3,T}) -> U

The value `geometry` takes at `pos`, for a voxel of side lengths `voxel_size`.
"""
@inline sample(geometry::Geometry, pos::NTuple{3}, voxel_size::NTuple{3}) =
    sample(geometry.shapes, pos, voxel_size, geometry.background, geometry.aa)

# Recurses on the tuple type to stay type stable and inline on the GPU.
# A shape that evaluates to background falls through to the next one.
@inline sample(::Tuple{}, pos, voxel_size, background, antialiasing) = background
@inline function sample(shapes::Tuple, pos, voxel_size, background, antialiasing)
    val = aa(first(shapes), pos, voxel_size, background, antialiasing)
    return val !== background ? val : sample(Base.tail(shapes), pos, voxel_size, background, antialiasing)
end
