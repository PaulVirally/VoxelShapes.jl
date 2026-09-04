"""
    Geometry{U, S, A}

Immutable description of what to draw: a tuple of shapes, a background value,
and an anti-aliasing strategy.

A `Geometry` says nothing about where the samples are taken. Hand it to
[`rasterize`](@ref) together with a [`Region`](@ref) or a
[`CompositeGrid`](@ref) to get an array, or call [`sample`](@ref) at a single
position.

Shapes are stored as a heterogeneous `Tuple` rather than a `Vector` so that
`sample` compiles to a type stable, GPU safe fold over them. They are evaluated
in order and the first one whose containment test succeeds claims the position.

The background value doubles as a transparency sentinel. A shape that covers a
position but evaluates to a value equal to `background` is treated as
transparent there, and evaluation falls through to the next shape. Filling with
the background value therefore punches a hole rather than painting that value
over a lower layer.

Use [`add_shape`](@ref) to append a shape, which returns a new `Geometry`.

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

A vector is converted to a tuple, so the resulting `Geometry` is `isbits` whenever
its shapes are.
"""
Geometry(shapes::Tuple, background) = Geometry(shapes, background, NoAntiAliasing())
Geometry(shapes::AbstractVector, background, antialiasing::AbstractAntiAliasing=NoAntiAliasing()) =
    Geometry(Tuple(shapes), background, antialiasing)

Base.eltype(geometry::Geometry) = typeof(geometry.background)

"""
    add_shape(geometry, shape) -> Geometry

Return a new `Geometry` with `shape` appended to the end of the shape tuple.

`Geometry` is immutable and its shapes live in a `Tuple`, so this is the functional
replacement for `push!`.
"""
add_shape(geometry::Geometry, shape) = Geometry((geometry.shapes..., shape), geometry.background, geometry.aa)

"""
    sample(geometry, pos::NTuple{3,T}, voxel_size::NTuple{3,T}) -> U

The value `geometry` takes at `pos`, for a voxel of side lengths `voxel_size`.

This is the same computation [`rasterize`](@ref) performs at every voxel center.
Anti-aliasing needs the voxel size because it samples inside the voxel, so
`voxel_size` matters even though the result is a value at a point. It can be
passed to Gila as a function of position with `pos -> sample(geometry, pos, vs)`.
"""
@inline sample(geometry::Geometry, pos::NTuple{3}, voxel_size::NTuple{3}) =
    sample(geometry.shapes, pos, voxel_size, geometry.background, geometry.aa)

# The fold over the shape tuple recurses on the tuple type so that it stays type
# stable and inlines on the GPU. The background acts as the transparency
# sentinel, so a shape that evaluates to it falls through to the next one.
@inline sample(::Tuple{}, pos, voxel_size, background, antialiasing) = background
@inline function sample(shapes::Tuple, pos, voxel_size, background, antialiasing)
    val = aa(first(shapes), pos, voxel_size, background, antialiasing)
    return val !== background ? val : sample(Base.tail(shapes), pos, voxel_size, background, antialiasing)
end
