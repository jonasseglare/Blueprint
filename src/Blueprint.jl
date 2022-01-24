module Blueprint

using FunctionalCollections
using LinearAlgebra
using Setfield

### Planes

const PlaneKey = Symbol

## Represents the plane normal*X = offset
struct Plane{T}
    normal::Vector{T}
    offset::T
end

function flip(plane::Plane{T}) where {T}
    return Plane{T}(-plane.normal, -plane.offset)
end

function plane_at_pos(normal::Vector{T}, pos::Vector{T}) where {T}
    Plane(normal, dot(pos, normal))
end

function pos_in_plane(plane::Plane{T}) where {T}
    lambda = plane.offset/dot(plane.normal, plane.normal)
    return lambda*plane.normal
end

function evaluate(plane::Plane{T}, pos::Vector{T}) where {T}
    return dot(plane.normal, pos) - plane.offset
end

function inside_halfspace(plane::Plane{T}, pos::Vector{T}) where {T}
    return evaluate(plane, pos) >= 0.0
end

function scale(x::T, plane::Plane{T}) where {T}
    return Plane{T}(x*plane.normal, x*plane.offset)
end

function normalize_plane(plane::Plane{T}) where {T}
    len = norm(plane.normal)
    return scale(1.0/len, plane)
end

## Line

struct ParameterizedLine{T}
    dir::Vector{T}
    pos::Vector{T}
end

function evaluate(line::ParameterizedLine{T}, lambda::T) where {T}
    return line.pos + lambda*line.dir
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




### Rigid transforms
struct RigidTransform{T}
    rotation::Matrix{T}
    translation::Vector{T}
end

function rigid_transform_from_xy_rotation(angle::Float64, dims::Integer)
    # #rotation = zeros(dims, dims)
    # translation = zeros(dims)
    cosa = cos(angle)
    sina = sin(angle)
    # return RigidTransform(rotation, translation)
    dst = identity_rigid_transform(dims)
    dst.rotation[1, 1] = cosa
    dst.rotation[2, 2] = cosa
    dst.rotation[1, 2] = -sina
    dst.rotation[2, 1] = sina
    return dst
end

function rigid_transform_from_rotation(R::Matrix{T}) where {T}
    (rows, cols) = size(R)
    @assert rows == cols
    return RigidTransform{T}(R, zeros(T, rows))
end

function rigid_transform_from_translation(v::Vector{Float64})
    n = length(v)
    return RigidTransform(Matrix(1.0I, n, n), v)
end

function transform_position(transform::RigidTransform{T}, position::Vector{T}) where {T}
    return transform.rotation*position + transform.translation
end

function transform_direction(transform::RigidTransform{T}, direction::Vector{T}) where {T}
    return transform.rotation*direction
end

function transform(transform::RigidTransform{T}, plane::Plane{T}) where {T}
    normal = transform_direction(transform, plane.normal)
    pos = transform_position(transform, pos_in_plane(plane))
    return plane_at_pos(normal, pos)
end

function identity_rigid_transform(N::Integer)
    RigidTransform(Matrix(1.0I, N, N), zeros(N))
end

function invert(transform::RigidTransform{T}) where {T}
    Rt = transform.rotation'
    return RigidTransform{T}(Rt, -Rt*transform.translation)
end

function compose(a::RigidTransform{T}, b::RigidTransform{T}) where {T}
    return RigidTransform{T}(a.rotation*b.rotation, a.rotation*b.translation + a.translation)
end

const identity_transform_3d = identity_rigid_transform(3)

function ordered_pair(a::PlaneKey, b::PlaneKey)::Tuple{PlaneKey,PlaneKey}
    if a < b
        (a, b)
    else
        (b, a)
    end
end

function ordered_triplet(a::PlaneKey, b::PlaneKey, c::PlaneKey)
    Tuple(sort([a, b, c]))
end

struct PlaneBound
    value::Float64
    key::PlaneKey
end

function update_plane_bound(
    bound::Union{PlaneBound,Nothing}, key::PlaneKey, value::Float64, cmp)
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


function update_line_bounds(dst::LineBounds, key::PlaneKey, plane::Plane{Float64})::LineBounds
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

function compute_bounded_lines(plane_map::PersistentHashMap{PlaneKey,Plane{Float64}})
    # List all lines between pairs of intersecting planes
    all_plane_intersections = Dict{Tuple{PlaneKey,PlaneKey}, ParameterizedLine{Float64}}()
    for x in plane_map
        for y in plane_map
            line = intersect(x.second, y.second)
            if line != nothing
                all_plane_intersections[ordered_pair(x.first, y.first)] = line
            end
        end
    end

    plane_intersections = Dict{Tuple{PlaneKey,PlaneKey}, LineBounds}()
    
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
    
    return phmap(collect(plane_intersections))
end

function shadowed_by(a::Plane{T}, b::Plane{T}) where {T}
    return a.normal == b.normal && a.offset < b.offset
end

function shadowed_by(a::Plane{T}, plane_map::AbstractDict{PlaneKey, Plane{T}}) where {T}
    for (sym, plane) in plane_map
        if shadowed_by(a, plane)
            return true
        end
    end
    return false
end

function remove_shadowed_planes(plane_map::AbstractDict{PlaneKey,Plane{Float64}})
    dst = Dict{PlaneKey,Plane{Float64}}()
    for (sym, plane) in plane_map
        if !shadowed_by(plane, plane_map)
            dst[sym] = plane
        end
    end
    return phmap(collect(dst))
end

# A polyhedron is defined by the half spaces of a set of planes. The half space of a plane is the space in
# the positive direction of the normal.
struct Polyhedron
    planes::PersistentHashMap{PlaneKey,Plane{Float64}}
    bounded_lines::PersistentHashMap{Tuple{PlaneKey,PlaneKey}, LineBounds}
    corners::Dict{Tuple{PlaneKey,PlaneKey,PlaneKey}, Vector{Float64}}
end

function visit_bounded_line_corner!(
    a::PlaneKey, b::PlaneKey, bds::LineBounds,
    pb::Union{PlaneBound, Nothing},
    dst::Dict{Tuple{PlaneKey, PlaneKey, PlaneKey}, Vector{Float64}})
    if pb != nothing
        dst[ordered_triplet(a, b, pb.key)] = evaluate(bds.line, pb.value)
    end
end

function collect_corners(
    bounded_lines::PersistentHashMap{Tuple{PlaneKey,PlaneKey}, LineBounds})
    dst = Dict{Tuple{PlaneKey, PlaneKey, PlaneKey}, Vector{Float64}}()
    for ((a, b), bds) in bounded_lines
        visit_bounded_line_corner!(a, b, bds, bds.lower, dst)
        visit_bounded_line_corner!(a, b, bds, bds.upper, dst)
    end
    return dst #phmap(collect(dst))
end

function polyhedron_from_planes(plane_map::AbstractDict{PlaneKey,Plane{Float64}})::Polyhedron
    plane_map = remove_shadowed_planes(plane_map)
    bounded_lines = compute_bounded_lines(plane_map)
    corners = collect_corners(bounded_lines)
    return Polyhedron(plane_map, bounded_lines, corners)
end

function add_planes(polyhedron::Polyhedron, plane_map::AbstractDict{PlaneKey,Plane{Float64}})::Polyhedron
    return polyhedron_from_planes(merge(polyhedron.planes, plane_map))
end

function add_plane(polyhedron::Polyhedron, key::PlaneKey, value::Plane{Float64})::Polyhedron
    return add_planes(polyhedron, Dict(key => value))
end

function intersect(a::Polyhedron, b::Polyhedron)
    return add_planes(a, b.planes)
end

function transform(rigid_transform::RigidTransform{Float64}, polyhedron::Polyhedron)
    dst = Dict{PlaneKey, Plane{Float64}}()
    for (k, v) in polyhedron.planes
        dst[k] = transform(rigid_transform, v)
    end
    return polyhedron_from_planes(dst)
end

function bounding_points(polyhedron::Polyhedron)::AbstractVector{Vector{Float64}}
    return [map(bdl -> bdl.line.pos, values(polyhedron.bounded_lines)); collect(values(polyhedron.corners))]
end

### Beams

# X-width 
# Y-height # Usually, X < Y
# Z-length

const world_up = [0.0, 0.0, 1.0]
const local_beam_dir = [0.0, 0.0, 1.0]
const beam_X_lower = :beam_X_lower
const beam_X_upper = :beam_X_upper
const beam_Y_lower = :beam_Y_lower
const beam_Y_upper = :beam_Y_upper

struct BeamSpecs
    Xsize::Real
    Ysize::Real

    # The Z direction is the direction of the beam
end

struct Beam
    transform::RigidTransform{Float64}
    specs::BeamSpecs
    polyhedron::Polyhedron
end

function plane_at_dim(dim::Integer, pos::Float64)
    normal = [0.0, 0.0, 0.0]
    normal[dim] = 1.0
    return Plane{Float64}(normal, pos)
end

function new_beam(beam_specs::BeamSpecs)
    return Beam(identity_transform_3d,
    beam_specs,
    polyhedron_from_planes(Dict(beam_X_lower => plane_at_dim(1, 0.0),
                                beam_X_upper => flip(plane_at_dim(1, beam_specs.Xsize)),
                                beam_Y_lower => plane_at_dim(2, 0.0),
                                beam_Y_upper => flip(plane_at_dim(2, beam_specs.Ysize)))))
end

function orthogonalize(target, reference)
    refhat = normalize(reference)
    return target - dot(target, refhat)*refhat
end

function orthonormalize(target, reference)
    return normalize(orthogonalize(target, reference))
end

function orient_beam_transform(world_dir::Vector{Float64}, up_in_beam_coordinates::Vector{Float64})
    dir = normalize(world_dir)
    local_ortho_up = orthonormalize(up_in_beam_coordinates, local_beam_dir)
    local_complement = cross(local_ortho_up, local_beam_dir)
    world_ortho_up = orthonormalize(world_up, dir)
    world_complement = cross(world_ortho_up, dir)

    A = [local_ortho_up local_complement local_beam_dir]
    B = [world_ortho_up world_complement dir]
    T = rigid_transform_from_rotation(B*A')
    @assert 0 < det(T.rotation)
    return T
end

function transform(rigid_transform::RigidTransform{Float64}, beam::Beam)
    return Beam(compose(rigid_transform, beam.transform),
    beam.specs,
    transform(rigid_transform, beam.polyhedron))
end

function set_transform(beam::Beam, new_transform::RigidTransform{Float64})
    old_transform = beam.transform
    transform_to_apply = compose(new_transform, invert(old_transform))
    return transform(transform_to_apply, beam)
end

function orient_beam(beam::Beam, world_dir::Vector{Float64}, up_in_beam_coordinates::Vector{Float64})
    return set_transform(beam, orient_beam_transform(world_dir, up_in_beam_coordinates))
end

function bounding_points(beam::Beam)::AbstractVector{Vector{Float64}}
    return bounding_points(beam.polyhedron)
end

function push_against(plane::Plane{Float64}, beam::Beam)
    plane = normalize_plane(plane)
    pts = bounding_points(beam)
    if length(pts) == 0
        return beam
    end
    shifts = map(pt -> evaluate(plane, pt), pts)
    shift = minimum(shifts)
    translation = rigid_transform_from_translation(-shift*plane.normal)
    return transform(translation, beam)
end

struct NamedPlane
    name::PlaneKey
    plane::Plane{Float64}
end

function cut(plane::NamedPlane, beam::Beam)
    return @set beam.polyhedron = add_plane(beam.polyhedron, plane.name, plane.plane)
end

# What should be included with `using`.
export plane_at_pos, BeamSpecs, beam_factory, new_beam

end # module
