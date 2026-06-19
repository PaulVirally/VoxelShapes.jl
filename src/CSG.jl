module CSG

export UnionShape, IntersectionShape, DifferenceShape, ComplementShape
export csg_union, csg_intersect, csg_diff, csg_complement

using ..Types

"""
    UnionShape{A, B} <: AbstractFillableShape

CSG union of two shapes. A point is inside if it is inside `a` or `b`.

The SDF is `min(sdf(a, p), sdf(b, p))`. Fill delegates to `a` when the point
is inside `a`, otherwise to `b`. The interpolation scheme of `a` is used.
`has_exact_sdf` is `true` only when both operands have exact SDFs.
"""
struct UnionShape{A<:AbstractFillableShape, B<:AbstractFillableShape} <: AbstractFillableShape
    a::A
    b::B
end

"""
    IntersectionShape{A, B} <: AbstractFillableShape

CSG intersection of two shapes. A point is inside if it is inside both `a` and `b`.

The SDF is `max(sdf(a, p), sdf(b, p))`. Fill always delegates to `a`.
`has_exact_sdf` is `true` only when both operands have exact SDFs.
"""
struct IntersectionShape{A<:AbstractFillableShape, B<:AbstractFillableShape} <: AbstractFillableShape
    a::A
    b::B
end

"""
    DifferenceShape{A, B} <: AbstractFillableShape

CSG difference `a \\ b`. A point is inside if it is inside `a` but not `b`.

The SDF is `max(sdf(a, p), -sdf(b, p))`. Fill always delegates to `a`.
`has_exact_sdf` is `true` only when both operands have exact SDFs.
"""
struct DifferenceShape{A<:AbstractFillableShape, B<:AbstractFillableShape} <: AbstractFillableShape
    a::A
    b::B
end

"""
    ComplementShape{A} <: AbstractFillableShape

Complement of a shape. A point is inside if it is outside `a`.

The SDF is `-sdf(a, p)`. Fill delegates to `a` (intended for use inside
an [`IntersectionShape`](@ref), not direct rendering).
"""
struct ComplementShape{A<:AbstractFillableShape} <: AbstractFillableShape
    a::A
end

"""
    csg_union(a, b) -> UnionShape

Construct the CSG union of two shapes. See [`UnionShape`](@ref).
"""
csg_union(a, b) = UnionShape(a, b)

"""
    csg_intersect(a, b) -> IntersectionShape

Construct the CSG intersection of two shapes. See [`IntersectionShape`](@ref).
"""
csg_intersect(a, b) = IntersectionShape(a, b)

"""
    csg_diff(a, b) -> DifferenceShape

Construct the CSG difference `a \\ b`. See [`DifferenceShape`](@ref).
"""
csg_diff(a, b) = DifferenceShape(a, b)

"""
    csg_complement(a) -> ComplementShape

Construct the complement of a shape. See [`ComplementShape`](@ref).
"""
csg_complement(a) = ComplementShape(a)

# interpolation: delegate to the primary operand (a)
Types.interpolation(s::UnionShape) = interpolation(s.a)
Types.interpolation(s::IntersectionShape) = interpolation(s.a)
Types.interpolation(s::DifferenceShape) = interpolation(s.a)
Types.interpolation(s::ComplementShape) = interpolation(s.a)

# Base.in
Base.in(p, s::UnionShape) = p in s.a || p in s.b
Base.in(p, s::IntersectionShape) = p in s.a && p in s.b
Base.in(p, s::DifferenceShape) = p in s.a && !(p in s.b)
Base.in(p, s::ComplementShape) = !(p in s.a)

# Base.fill: delegate to primary operand (a)
Base.fill(s::UnionShape, vc, vs) = Tuple(vc) in s.a ? fill(s.a, vc, vs) : fill(s.b, vc, vs)
Base.fill(s::IntersectionShape, vc, vs) = fill(s.a, vc, vs)
Base.fill(s::DifferenceShape, vc, vs) = fill(s.a, vc, vs)
# ComplementShape is meant to be used inside IntersectionShape, not rendered directly.
Base.fill(s::ComplementShape, vc, vs) = fill(s.a, vc, vs)

# SDFs: classic CSG distance combinators
Types.sdf(s::UnionShape, p) = min(sdf(s.a, p), sdf(s.b, p))
Types.sdf(s::IntersectionShape, p) = max(sdf(s.a, p), sdf(s.b, p))
Types.sdf(s::DifferenceShape, p) = max(sdf(s.a, p), -sdf(s.b, p))
Types.sdf(s::ComplementShape, p) = -sdf(s.a, p)

# has_exact_sdf: (conservative) only true when both operands have exact SDFs
Types.has_exact_sdf(s::UnionShape) = has_exact_sdf(s.a) && has_exact_sdf(s.b)
Types.has_exact_sdf(s::IntersectionShape) = has_exact_sdf(s.a) && has_exact_sdf(s.b)
Types.has_exact_sdf(s::DifferenceShape) = has_exact_sdf(s.a) && has_exact_sdf(s.b)
Types.has_exact_sdf(s::ComplementShape) = has_exact_sdf(s.a)

end # module
