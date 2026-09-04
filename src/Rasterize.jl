#= The voxel center of index idx is origin + (idx - 1/2) * voxel_size, the same
convention Gila uses for its cell centers. Only the converted origin and voxel
size reach the kernel: a Region holds exact rationals and is not isbits, so it
stays on the host. =#
@kernel function fill_voxel!(arr, geometry::Geometry, origin::NTuple{3}, voxel_size::NTuple{3})
    idx = @index(Global, Cartesian)
    half = one(eltype(voxel_size)) / 2
    pos = ntuple(dir -> origin[dir] + (idx[dir] - half) * voxel_size[dir], 3)
    arr[idx] = sample(geometry, pos, voxel_size)
end

"""
    rasterize!(arr, geometry, region::Region; precision=Float64) -> arr

Fill the three dimensional array `arr` with the values of `geometry` at the voxel
centers of `region`.

The backend is taken from the type of `arr`, so an `Array` runs on the CPU and a
`CuArray` on the GPU. `arr` must have `size(region)` elements.

# Arguments
- `arr`: the destination array, of size `size(region)`
- `geometry`: the [`Geometry`](@ref) to sample
- `region`: the [`Region`](@ref) giving the voxel centers
- `precision::Type{<:AbstractFloat}=Float64`: the floating point type the
  region's exact rational geometry is converted to before it reaches the kernel.
  It has to match the `T` of the shapes in `geometry`, so `Float32` shapes need
  `precision=Float32`.
"""
function rasterize!(arr::AbstractArray{U,3}, geometry::Geometry, region::Region;
                    precision::Type{<:AbstractFloat}=Float64) where {U}
    if size(arr) != size(region)
        throw(DimensionMismatch("The destination array is $(size(arr)) but the region has $(size(region)) voxels."))
    end
    kernel = fill_voxel!(KernelAbstractions.get_backend(arr))
    kernel(arr, geometry, precision.(lower_corner(region)), precision.(scale(region));
           ndrange=size(arr))
    return arr
end

"""
    rasterize(geometry, region::Region; precision=Float64) -> Array{U,3}
    rasterize(geometry, region::Region, ::Type{CuArray}; precision=Float64) -> CuArray{U,3}
    rasterize(geometry, grid::CompositeGrid; precision=Float64) -> CompositeField
    rasterize(geometry, grid::CompositeGrid, ::Type{CuArray}; precision=Float64) -> CompositeField

Sample `geometry` at every voxel center of a grid.

Over a [`Region`](@ref) the result is a dense array of size `size(region)`.
Over a [`CompositeGrid`](@ref) it is a [`CompositeField`](@ref), which holds one
flat vector in the region layout order of the grid. Passing `CuArray` as the
third argument puts the result on the GPU.

# Arguments
- `geometry`: the [`Geometry`](@ref) to sample
- `region` or `grid`: where to sample
- `precision::Type{<:AbstractFloat}=Float64`: the floating point type the exact
  rational geometry is converted to, which must match the `T` of the shapes
"""
function rasterize(geometry::Geometry, region::Region; precision::Type{<:AbstractFloat}=Float64)
    arr = Array{eltype(geometry),3}(undef, size(region))
    return rasterize!(arr, geometry, region; precision=precision)
end

function rasterize(geometry::Geometry, region::Region, ::Type{CuArray};
                   precision::Type{<:AbstractFloat}=Float64)
    arr = CuArray{eltype(geometry),3}(undef, size(region))
    return rasterize!(arr, geometry, region; precision=precision)
end

"""
    CompositeField{U, V}

A scalar field over a [`CompositeGrid`](@ref), stored as one flat vector.

The regions of a composite grid have different voxel counts, so a field over
them cannot be a single dense array. Instead the values of region `i` occupy the
contiguous block `offsets[i]+1 : offsets[i+1]` of `data`, column major inside the
block, with the regions in the order [`regions`](@ref) reports them. That is
exactly the layout Gila expects, so `vec(field)` can be handed to
`GlaCmpVol`-based operators directly, and `collect(eachregion(field))` gives the
per-region tensor form.

Use [`regionview`](@ref) for the values of one region as a three dimensional
array, and [`regrid`](@ref) for a single uniform array over the whole grid.

# Fields
- `data::V`: the values, flat, one entry per voxel in region layout order
- `grid::CompositeGrid`: the grid the field lives on
- `offsets::Vector{Int}`: block boundaries, `nregions(grid) + 1` entries starting at 0
"""
struct CompositeField{U, V<:AbstractVector{U}}
    data::V
    grid::CompositeGrid
    offsets::Vector{Int}
end

"""
    CompositeField(data, grid)

Wrap a flat vector of values as a field over `grid`, computing the region block
offsets from the grid.
"""
function CompositeField(data::AbstractVector{U}, grid::CompositeGrid) where {U}
    if length(data) != length(grid)
        throw(DimensionMismatch("The data vector has $(length(data)) entries but the grid has $(length(grid)) voxels."))
    end
    offsets = zeros(Int, nregions(grid) + 1)
    for (idx, reg) in enumerate(regions(grid))
        offsets[idx + 1] = offsets[idx] + length(reg)
    end
    return CompositeField{U,typeof(data)}(data, grid, offsets)
end

function rasterize(geometry::Geometry, grid::CompositeGrid; precision::Type{<:AbstractFloat}=Float64)
    field = CompositeField(Vector{eltype(geometry)}(undef, length(grid)), grid)
    return rasterize!(field, geometry; precision=precision)
end

function rasterize(geometry::Geometry, grid::CompositeGrid, ::Type{CuArray};
                   precision::Type{<:AbstractFloat}=Float64)
    field = CompositeField(CuVector{eltype(geometry)}(undef, length(grid)), grid)
    return rasterize!(field, geometry; precision=precision)
end

"""
    rasterize!(field::CompositeField, geometry; precision=Float64) -> field

Fill an existing [`CompositeField`](@ref) with the values of `geometry`, one kernel
launch per region.
"""
function rasterize!(field::CompositeField, geometry::Geometry;
                    precision::Type{<:AbstractFloat}=Float64)
    for idx in 1:nregions(field)
        rasterize!(regionview(field, idx), geometry, regions(field)[idx]; precision=precision)
    end
    return field
end

Base.vec(field::CompositeField) = field.data
Base.eltype(::CompositeField{U}) where {U} = U
Base.length(field::CompositeField) = length(field.data)

Grids.regions(field::CompositeField) = regions(field.grid)
Grids.nregions(field::CompositeField) = nregions(field.grid)

"""
    regionview(field::CompositeField, idx::Integer)

The values of region `idx` of `field`, reshaped to `cells(regions(field)[idx])`.

The result is a view into `vec(field)`, so writing to it writes through to the
field.
"""
function regionview(field::CompositeField, idx::Integer)
    block = (field.offsets[idx] + 1):field.offsets[idx + 1]
    return reshape(view(field.data, block), cells(field.grid.regions[idx]))
end

"""
    eachregion(field::CompositeField)

Iterate over the regions of `field` as three dimensional views, in the order
[`regions`](@ref) reports them.
"""
eachregion(field::CompositeField) = (regionview(field, idx) for idx in 1:nregions(field))

"""
    regrid(field::CompositeField, scale=finest(field.grid)) -> Array{U,3}

Resample `field` onto one uniform array covering the bounding box of its grid.

Every voxel of the result takes the value of the source voxel that contains its
center, so a coarse region is replicated over the finer output voxels it covers.
This is the cheapest way to look at a composite field, but it costs the memory a
uniform grid at that resolution needs.

`regrid` runs on the CPU. A field on the GPU is copied to the host first.

# Arguments
- `field`: the field to resample
- `scale=finest(field.grid)`: the voxel side lengths of the result

# Returns
- An `Array{U,3}` spanning `bounding_box(field.grid)`

# Throws
- `ArgumentError`: if `scale` does not divide the voxel size of every region, or
  does not divide the bounding box into a whole number of voxels
"""
function regrid(field::CompositeField, scl::NTuple{3,Union{Integer,Rational}}=finest(field.grid))
    scl = Rational{Int}.(scl)
    grid = field.grid
    for (idx, reg) in enumerate(regions(grid))
        ratio = scale(reg) .// scl
        if any(.!isinteger.(ratio))
            bad = findall(collect(.!isinteger.(ratio)))
            throw(ArgumentError("Region $idx has voxel size $(scale(reg)), which is not a whole number of $(scl) voxels in dimension(s) $(bad). Pick a scale that divides every region scale, such as finest(grid) = $(finest(grid))."))
        end
    end
    lwr, upr = bounding_box(grid)
    span = (upr .- lwr) .// scl
    if any(.!isinteger.(span))
        bad = findall(collect(.!isinteger.(span)))
        throw(ArgumentError("The grid spans $(upr .- lwr), which is not a whole number of $(scl) voxels in dimension(s) $(bad)."))
    end
    num = ntuple(dir -> Int(span[dir]), 3)

    data = Array(field.data)
    out = Array{eltype(field),3}(undef, num)
    regs = regions(grid)
    corners = [lower_corner(reg) for reg in regs]
    tops = [upper_corner(reg) for reg in regs]
    for idx in CartesianIndices(out)
        pos = ntuple(dir -> lwr[dir] + (idx[dir] - 1 // 2) * scl[dir], 3)
        # The regions tile the box, so exactly one of them contains this center.
        # Cell centers never land on a region boundary because scl divides every
        # region scale, hence the half open test.
        for (reg_idx, reg) in enumerate(regs)
            if all(corners[reg_idx] .<= pos) && all(pos .< tops[reg_idx])
                cel = ntuple(dir -> Int(fld(pos[dir] - corners[reg_idx][dir], scale(reg)[dir])) + 1, 3)
                flat = LinearIndices(cells(reg))[cel...]
                out[idx] = data[field.offsets[reg_idx] + flat]
                break
            end
        end
    end
    return out
end

Base.:(==)(fldA::CompositeField, fldB::CompositeField) =
    fldA.grid == fldB.grid && fldA.data == fldB.data

function Base.show(io::IO, field::CompositeField)
    num = nregions(field)
    print(io, "CompositeField{", eltype(field), "} over ", num, " region",
          num == 1 ? "" : "s", ", ", length(field), " voxels")
end
Base.show(io::IO, ::MIME"text/plain", field::CompositeField) = show(io, field)
