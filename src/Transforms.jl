module Transforms

using StaticArrays
using LinearAlgebra: norm

export Rotated

using ..Types

"""
    Rotated{S, T} <: AbstractFillableShape

Wrapper that rotates any `AbstractFillableShape` about a pivot point.

# Fields
- `inner`: the wrapped shape
- `R`: world-to-local rotation matrix, applied as `R * (p - pivot) + pivot`
- `pivot`: world-space point rotated about
"""
struct Rotated{S<:AbstractFillableShape, T} <: AbstractFillableShape
    inner::S
    R::SMatrix{3,3,T,9}
    pivot::SVector{3,T}
end

# Constructor from explicit SMatrix
function Rotated(inner::S, R::SMatrix{3,3,T,9}) where {S<:AbstractFillableShape, T}
    piv = SVector{3,T}(center(inner))
    Rotated{S,T}(inner, R, piv)
end

function Rotated(inner::S, R::SMatrix{3,3,T,9}, pivot::NTuple{3,T}) where {S<:AbstractFillableShape, T}
    Rotated{S,T}(inner, R, SVector{3,T}(pivot))
end

# Constructor from Euler angles (intrinsic ZYX / Tait-Bryan)
function Rotated(inner::S, angles::NTuple{3,T}) where {S<:AbstractFillableShape, T}
    αx, αy, αz = angles
    Rx = SMatrix{3,3,T,9}(1, 0, 0,
                           0, cos(αx), sin(αx),
                           0, -sin(αx), cos(αx))
    Ry = SMatrix{3,3,T,9}(cos(αy), 0, -sin(αy),
                           0, 1, 0,
                           sin(αy), 0, cos(αy))
    Rz = SMatrix{3,3,T,9}(cos(αz), sin(αz), 0,
                           -sin(αz), cos(αz), 0,
                           0, 0, 1)
    R = Rx * Ry * Rz
    Rotated(inner, R)
end

# Constructor from axis-angle
function Rotated(inner::S, axis::NTuple{3,T}, angle::T) where {S<:AbstractFillableShape, T}
    u = SVector{3,T}(axis) / norm(SVector{3,T}(axis))
    c, s = cos(angle), sin(angle)
    ux, uy, uz = u[1], u[2], u[3]
    R = SMatrix{3,3,T,9}(
        c + ux^2*(1-c),    ux*uy*(1-c) + uz*s, ux*uz*(1-c) - uy*s,
        uy*ux*(1-c) - uz*s, c + uy^2*(1-c),    uy*uz*(1-c) + ux*s,
        uz*ux*(1-c) + uy*s, uz*uy*(1-c) - ux*s, c + uz^2*(1-c)
    )
    Rotated(inner, R)
end

@inline _to_local(p::NTuple{3,T}, w::Rotated{S,T}) where {S,T} =
    Tuple(w.R * (SVector{3,T}(p) - w.pivot) + w.pivot)

Types.interpolation(w::Rotated) = interpolation(w.inner)

Base.in(point::NTuple{3,T}, w::Rotated{S,T}) where {S,T} =
    _to_local(point, w) in w.inner

Base.fill(w::Rotated{S,T}, vc::NTuple{3,T}, vs::NTuple{3,T}) where {S,T} =
    fill(w.inner, _to_local(vc, w), vs)

# Rotation is isometric, so SDF distances are preserved exactly.
Types.sdf(w::Rotated{S,T}, point::NTuple{3,T}) where {S,T} =
    sdf(w.inner, _to_local(point, w))

Types.has_exact_sdf(w::Rotated) = has_exact_sdf(w.inner)

"""
    bounding_box(shape::Rotated) -> (lower, upper)

All-infinite if the inner box is infinite along any axis. Otherwise it's the
componentwise min/max of the inner box's 8 corners mapped to world space
(`world = R' * (c - pivot) + pivot`).
"""
function Types.bounding_box(w::Rotated{S,T}) where {S,T}
    lo, hi = bounding_box(w.inner)
    if any(isinf, lo) || any(isinf, hi)
        inf = T(Inf)
        return ((-inf, -inf, -inf), (inf, inf, inf))
    end
    lo3 = SVector{3,T}(lo[1], lo[2], lo[3])
    hi3 = SVector{3,T}(hi[1], hi[2], hi[3])
    corners = ntuple(8) do k
        x = (k - 1) & 1 == 0 ? lo3[1] : hi3[1]
        y = (k - 1) & 2 == 0 ? lo3[2] : hi3[2]
        z = (k - 1) & 4 == 0 ? lo3[3] : hi3[3]
        w.R' * (SVector{3,T}(x, y, z) - w.pivot) + w.pivot
    end
    lower = ntuple(i -> minimum(c[i] for c in corners), 3)
    upper = ntuple(i -> maximum(c[i] for c in corners), 3)
    return (lower, upper)
end

end # module
