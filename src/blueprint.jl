module blueprint

using FunctionalCollections
using LinearAlgebra
using Setfield

### Planes

## Represents the plane normal*X = offset
struct Plane{T}
    normal::Vector{T}
    offset::T
end

function plane_at_pos(normal::Vector{T}, pos::Vector{T}) where {T}
    Plane(normal, dot(pos, normal))
end

function evaluate(plane::Plane{T}, pos::Vector{T}) where {T}
    return dot(plane.normal, pos) - plane.offset
end

function inside_halfspace(plane::Plane{T}, pos::Vector{T}) where {T}
    return evaluate(plane, pos) >= 0.0
end    

## Line

struct ParameterizedLine{T}
    dir::Vector{T}
    pos::Vector{T}
end

function intersect(a::Plane{T}, b::Plane{T})::Union{ParameterizedLine{T}, Nothing} where {T}
    dir = cross(a.normal, b.normal)
    if norm(dir) <= 0.0
        return nothing
    end

    # Find the point that satisfies the plane equations and 
    origin = (hcat(a.normal, b.normal, dir)')\[a.offset, b.offset, 0.0]
    return ParameterizedLine{T}(dir, origin)
end

struct PlaneLineIntersection{T}
    normal_dir_dot::T
    lambda::T
end

function exists(intersection::PlaneLineIntersection{T})::Bool where {T}
    return intersection.normal_dir_dot != 0.0
end

function intersect(a::Plane{T}, b::ParameterizedLine{T})::PlaneLineIntersection where {T}
    denom = dot(a.normal, b.dir)
    return PlaneLineIntersection{T}(denom, (a.offset - dot(a.normal, b.pos))/denom)
end



### Beams

# X-width 
# Y-height # Usually, X < Y
# Z-length

struct BeamSpecs
    width::Real # along X
    height::Real # along Y
end

function conventional_specs(specs::BeamSpecs)
    return specs.width <= specs.height
end

### Rigid transforms
struct RigidTransform{T}
    rotation::Matrix{T}
    translation::Vector{T}
end

function identity_rigid_transform(N::Integer)
    RigidTransform(Matrix(1.0I, N, N), zeros(N))
end

const identity_3d_transform = identity_rigid_transform(3)
    
### Beams

mutable struct BeamFactory
    prefix::String
    specs::BeamSpecs
    counter::Integer
end

function beam_factory(prefix::String, specs::BeamSpecs)
    if !conventional_specs(specs)
        print("Specs " + prefix + " are not conventional")
    end
    BeamFactory(prefix, specs, 0)
end

# A polyhedron is defined by the half spaces of a set of planes. The half space of a plane is the space in
# the positive direction of the normal.
struct Polyhedron
end

function ordered_pair(a::Symbol, b::Symbol)::Tuple{Symbol,Symbol}
    if a < b
        (a, b)
    else
        (b, a)
    end
end

struct PlaneBound
    value::Float64
    key::Symbol
end

function update_plane_bound(
    bound::Union{PlaneBound,Nothing}, key::Symbol, value::Float64, cmp)
    if bound == nothing
        return PlaneBound(value, key)
    elseif cmp(bound.value, value)
        return PlaneBound(value, key)
    else
        return bound
    end
end

struct LineBounds
    line::ParameterizedLine{Float64}
    exists::Bool

    # Nothing means unlimited
    lower::Union{PlaneBound,Nothing}
    upper::Union{PlaneBound,Nothing}
end

function initialize_line_bounds(line::ParameterizedLine{Float64})
    return LineBounds(line, true, nothing, nothing)
end

function check_existence(bds::LineBounds)
    if bds.lower != nothing && bds.upper != nothing && bds.lower.value >= bds.upper.value
        return @set bds.exists = false
    else
        return bds
    end
end

function or_nothing(a, b)
    if a == nothing
        return b
    else
        return a
    end
end


function update_line_bounds(dst::LineBounds, key::Symbol, plane::Plane{Float64})::LineBounds
    if !dst.exists
        return false
    end

    intersection = intersect(plane, dst.line)
    if exists(intersection)
        # If the line direction points in opposite direction as normal,
        # the plane is upper bounding the line.
        is_upper = intersection.normal_dir_dot < 0

        l = intersection.lambda
        return check_existence(
        if is_upper
            @set dst.upper = update_plane_bound(dst.upper, key, l, >)
        else
            @set dst.lower = update_plane_bound(dst.lower, key, l, <)
        end)
    else
        if inside_halfspace(plane, dst.line.pos)
            return dst
        else
            return @set dst.exists = false
        end
    end
end

struct PolyhedronSettings
    marg::Float64
end

function default_polyhedron_settings()
    return PolyhedronSettings(1.0e-6)
end

function compute_bounded_lines(plane_map::PersistentHashMap{Symbol,Plane{Float64}})
    # List all lines between pairs of intersecting planes
    all_plane_intersections = Dict{Tuple{Symbol,Symbol}, ParameterizedLine{Float64}}()
    for x in plane_map
        for y in plane_map
            line = intersect(x.second, y.second)
            if line != nothing
                all_plane_intersections[ordered_pair(x.first, y.first)] = line
            end
        end
    end

    plane_intersections = Dict{Tuple{Symbol,Symbol}, LineBounds}()
    
    # For all lines, remove lines that don't exist because of other planes
    for ((p0, p1), line) in all_plane_intersections
        is_included = true
        bounds = initialize_line_bounds(line)
        for (plane_key, plane) in plane_map
            if plane_key != p0 && plane_key != p1
                bounds = update_line_bounds(bounds, plane_key, plane)
            end
        end

        if bounds.exists
            plane_intersections[(p0, p1)] = bounds
        end
    end
    
    return plane_intersections
end

function remove_redundant_planes(plane_map::PersistentHashMap{Symbol,Plane{Float64}})
    return plane_map
end

function polyhedron_from_planes(plane_map::PersistentHashMap{Symbol,Plane{Float64}})
    plane_map = remove_redundant_planes(plane_map)
    bounded_lines = compute_bounded_lines(plane_map)
end

struct Beam
    name::String
    specs::BeamSpecs
    polyhedron::Polyhedron
end

function new_beam(name::String, specs::BeamSpecs)
    p = Polyhedron()
    Beam(name, specs, p)
end

function new_beam!(factory::BeamFactory)
    result = new_beam(string(factory.prefix, string(factory.counter)), factory.specs)
    factory.counter += 1
    result
end


# What should be included with `using`.
export plane_at_pos, BeamSpecs, beam_factory, new_beam

end # module
