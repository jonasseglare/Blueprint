module Blueprint

using FunctionalCollections
using LinearAlgebra
using Setfield

#using Luxor
import Luxor
const lx = Luxor

abstract type PhysicalObject end

### Planes

const PlaneKey = Symbol
const PlaneKeyTuple2 = Tuple{PlaneKey, PlaneKey}
const PlaneKeyTuple3 = Tuple{PlaneKey, PlaneKey, PlaneKey}

struct DefinedInterval{T}
    lower::T
    upper::T
end

function width(interval::DefinedInterval{T}) where {T}
    return interval.upper - interval.lower
end

const Interval = Union{DefinedInterval{T}, Nothing} where {T}

function undefined_interval() where {T}
    return nothing
end

function is_defined(x::Interval{T}) where {T}
    return x != nothing
end

function extend(interval::Interval{T}, x::T) where {T}
    if is_defined(interval)
        return DefinedInterval{T}(min(interval.lower, x), max(interval.upper, x))
    else
        return DefinedInterval{T}(x, x)
    end
end

struct BBox{T}
    intervals::Vector{Interval{T}}
end

function compute_bbox(points::AbstractVector{Vector{T}}) where {T}
    if length(points) == 0
        return BBox{T}([])
    end
    
    m = length(points)
    n = length(points[1])
    dst = Vector{Interval{T}}(undefined_interval(), n)
    for point in points
        for i in 1:n
            result = extend(dst[i], point[i])
            dst[i] = result
        end
    end
    return BBox{T}(dst)
end


## Represents the plane normal*X = offset
struct Plane{T}
    normal::Vector{T}
    offset::T
end

function flip(plane::Plane{T}) where {T}
    return Plane{T}(-plane.normal, -plane.offset)
end

function translate(plane::Plane{T}, amount::T) where {T}
    return Plane{T}(plane.normal, plane.offset + amount)
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

function project(plane::Plane{T}, X::Vector{T}) where {T}
    plane = normalize_plane(plane)
    return X - evaluate(plane, X)*plane.normal
end

function parallel_plane_distance(a::Plane{T}, b::Plane{T}) where {T}
    X = pos_in_plane(a)
    return norm(X - project(b, X))
end
    

function plane_at_dim(dim::Integer, pos::Float64)
    normal = [0.0, 0.0, 0.0]
    normal[dim] = 1.0
    return Plane{Float64}(normal, pos)
end

function facing_plane(plane::Plane{T}, distance) where {T}
    return flip(translate(plane, distance))
end


## Line

struct ParameterizedLine{T}
    dir::Vector{T}
    pos::Vector{T}
end

function evaluate(line::ParameterizedLine{T}, lambda::T) where {T}
    return line.pos + lambda*line.dir
end

function direct_like(line::ParameterizedLine{T}, dir::Vector{T}) where {T}
    if dot(line.dir, dir) < 0
        return ParameterizedLine{T}(-line.dir, line.pos)
    else
        return line
    end
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

function ordered_triplet(ab::Tuple{PlaneKey, PlaneKey}, c::PlaneKey)
    (a, b) = ab
    return ordered_triplet(a, b, c)
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
    planes::AbstractDict{PlaneKey,Plane{Float64}}
    bounded_lines::AbstractDict{Tuple{PlaneKey,PlaneKey}, LineBounds}
    corners::AbstractDict{Tuple{PlaneKey,PlaneKey,PlaneKey}, Vector{Float64}}
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

function remove_planes_without_corners(
    planes::AbstractDict{PlaneKey, Plane{Float64}},
    corners::Dict{Tuple{PlaneKey,PlaneKey,PlaneKey}, Vector{Float64}})

    if length(corners) == 0
        return planes
    end

    # If there exists at least one pair of planes that are not parallel, then
    # *all* planes that are part of the polyhedron will be have at least one corner.

    referred_planes = Set{PlaneKey}()
    for (k3, x) in corners
        for k in k3
            push!(referred_planes, k)
        end
    end
    
    dst = Dict{PlaneKey, Plane{Float64}}()
    for (k, p) in planes
        if k in referred_planes
            dst[k] = p
        end
    end
    return planes
end

function polyhedron_from_planes(plane_map::AbstractDict{PlaneKey,Plane{Float64}})::Polyhedron
    plane_map = remove_shadowed_planes(plane_map)
    bounded_lines = compute_bounded_lines(plane_map)
    corners = collect_corners(bounded_lines)
    return Polyhedron(remove_planes_without_corners(plane_map, corners), bounded_lines, corners)
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

function bounding_points_from_lines(polyhedron::Polyhedron)
    return map(bdl -> bdl.line.pos, values(polyhedron.bounded_lines))
end

function bounding_points_from_corners(polyhedron::Polyhedron)
    return collect(values(polyhedron.corners))
end

function bounding_points(polyhedron::Polyhedron)::AbstractVector{Vector{Float64}}
    corners = bounding_points_from_corners(polyhedron)
    if length(corners) == 0
        return bounding_points_from_lines(polyhedron)
    else
        return corners
    end
end

function mid_point(polyhedron::Polyhedron)
    pts = bounding_points(polyhedron)
    sum = pts[1]
    for pt in pts[2:end]
        sum += pt
    end
    return (1.0/length(pts))*sum
end

### Beams

# X-width (the shorter)
# Y-height (the longer) # Usually, X < Y
# Z-length

## Example:
#
#   Y
#   ^
#   |####
#   |####
#   |####
#   |####
#   |####
#   |####
#   o-----> X
#
const local_x_dir = [1.0, 0.0, 0.0]
const local_y_dir = [0.0, 1.0, 0.0]

const world_up = [0.0, 0.0, 1.0]
const local_beam_dir = [0.0, 0.0, 1.0]
const beam_X_lower = :beam_X_lower
const beam_X_upper = :beam_X_upper
const beam_Y_lower = :beam_Y_lower
const beam_Y_upper = :beam_Y_upper
const plane_keys_per_beam_dim = [(beam_X_lower, beam_X_upper), (beam_Y_lower, beam_Y_upper)]

struct RgbColor
    # in ranges 0..1
    red::Float64
    green::Float64
    blue::Float64
end

struct BeamSpecs
    Xsize::Real
    Ysize::Real

    color::RgbColor
    # The Z direction is the direction of the beam
end

default_beam_color = RgbColor(0.0, 0.0, 1.0)

function beam_specs(Xsize::Real, Ysize::Real)
    return BeamSpecs(Xsize, Ysize, default_beam_color)
end

function quadratic_beam_specs(size::Real)
    return beam_specs(size, size)
end

struct Annotation
    position::Vector{Float64}
    label::String
end


struct Beam <: PhysicalObject
    transform::RigidTransform{Float64}
    specs::BeamSpecs
    polyhedron::Polyhedron
    annotations::AbstractDict{PlaneKey, PersistentVector{Annotation}}
end

function new_beam(beam_specs::BeamSpecs)
    xlow = plane_at_dim(1, 0.0)
    ylow = plane_at_dim(2, 0.0)
    return Beam(identity_transform_3d,
    beam_specs,
    polyhedron_from_planes(Dict(beam_X_lower => xlow,
                                beam_X_upper => flip(plane_at_dim(1, beam_specs.Xsize)),
                                beam_Y_lower => ylow,
                                beam_Y_upper => flip(plane_at_dim(2, beam_specs.Ysize)))),
                Dict{PlaneKey, PersistentVector{Annotation}}())
end

function beam_dir(beam::Beam)
    return transform_direction(beam.transform, local_beam_dir)
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

function transform(rigid_transform::RigidTransform{Float64},
    src::AbstractDict{PlaneKey, PersistentVector{Annotation}})
    dst = Dict{PlaneKey, PersistentVector{Annotation}}()
    for (k, annotations) in src
        
        #dst[k] = pvec(map(annotation -> @set annotation.position = transform_position(rigid_transform, annotation.position), annotations))
        new_annots = [@set annotation.position = transform_position(
            rigid_transform, annotation.position) for annotation in annotations]
        dst[k] = pvec(new_annots)
    end
    return dst
end


function transform(rigid_transform::RigidTransform{Float64}, beam::Beam)
    return Beam(compose(rigid_transform, beam.transform),
    beam.specs,
    transform(rigid_transform, beam.polyhedron),
    transform(rigid_transform, beam.annotations))
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

function mid_point(beam::Beam)
    return mid_point(beam.polyhedron)
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

function contains_annotation(p::Polyhedron, k::PlaneKey, x::Vector{Float64})
    if !(haskey(p.planes, k))
        return false
    end
    
    for (k2, plane) in p.planes
        if k != k2
            if k != k2 && !(inside_halfspace(plane, x))
                return false
            end
        end
    end
    return true
end

function refilter_annotations(beam::Beam)
    dst = Dict{PlaneKey, PersistentVector{Annotation}}()
    for (plane_key, annotations) in beam.annotations
        filtered = filter(annotation -> contains_annotation(beam.polyhedron, plane_key, annotation.position), annotations)
        if 0 < length(filtered)
            dst[plane_key] = pvec(filtered)
        end
    end
    return @set beam.annotations = dst
end

function cut(plane::NamedPlane, beam::Beam)
    polyhedron = add_plane(beam.polyhedron, plane.name, plane.plane)
    return refilter_annotations(@set beam.polyhedron = polyhedron)
end

function compute_drilling_direction(first_beam::Beam, second_beam::Beam)
    first_dir = beam_dir(first_beam)
    second_dir = beam_dir(second_beam)
    dir = normalize(cross(first_dir, second_dir))
    a_pt = mid_point(first_beam)
    b_pt = mid_point(second_beam)
    diff = b_pt - a_pt
    if dot(dir, diff) >= 0
        return dir
    else
        return -dir
    end
end

function base_drilling_dim(beam::Beam, drilling_dir::Vector{Float64})
    local_drilling_dir = normalize(transform_direction(invert(beam.transform), drilling_dir))
    if abs(dot(local_drilling_dir, local_x_dir)) < abs(dot(local_drilling_dir, local_y_dir))
        return 1
    else
        return 2
    end
end

struct DrillingPlaneSpecs
    # How far from the bounding planes the drilling planes should be
    marg::Real
    
    # Number of planes to generate
    count::Integer 
end

function generate_drilling_planes(
    beam::Beam, specs::DrillingPlaneSpecs, drilling_dir::Vector{Float64})::Vector{Plane{Float64}}
    dim = base_drilling_dim(beam, drilling_dir)
    (lower_key, upper_key) = plane_keys_per_beam_dim[dim]
    lower_plane = beam.polyhedron.planes[lower_key]
    upper_plane = beam.polyhedron.planes[upper_key]
    dist = evaluate(lower_plane, pos_in_plane(upper_plane))
    lower_offset = specs.marg
    upper_offset = dist - specs.marg
    dst = Vector{Plane{Float64}}()
    if lower_offset >= upper_offset
        return dst
    end
    step = (upper_offset - lower_offset)/(specs.count - 1)
    return map(i -> translate(lower_plane, lower_offset + i*step), 0:(specs.count - 1))
end

struct Drill
    label::String
    line::ParameterizedLine{Float64}
end

function generate_drills(drilling_dir::Vector{Float64},
    a_planes::Vector{Plane{Float64}},
    b_planes::Vector{Plane{Float64}})::Vector{Drill}
    dst = Vector{Drill}()
    for a in a_planes
        for b in b_planes
            push!(dst, Drill("Drill", direct_like(intersect(a, b), drilling_dir)))
        end
    end
    return dst
end

function drill(beam::Beam, drills::AbstractVector{Drill})
    new_annotations = Dict{PlaneKey, PersistentVector{Annotation}}()
    for (plane_key, plane) in beam.polyhedron.planes
        plane = normalize_plane(plane)
        dst = get(beam.annotations, plane_key, pvec(Vector{Annotation}()))
        for drill in drills
            if dot(drill.line.dir, plane.normal) > 0
                intersection = intersect(plane, drill.line)
                if exists(intersection)
                    dst = push(dst, Annotation(evaluate(drill.line, intersection.lambda), drill.label))
                end
            end
        end
        new_annotations[plane_key] = dst
    end
    return refilter_annotations(@set beam.annotations = new_annotations)
end



function remove_from_tuple(x, k)
    return Tuple(filter(k2 -> k2 != k, x))
end

function plane_corner_keys(phd::Polyhedron, pk::PlaneKey)
    return [remove_from_tuple(k, pk) for (k, v) in phd.corners if pk in k]
end

function are_connected(a::Tuple{PlaneKey,PlaneKey}, b::Tuple{PlaneKey,PlaneKey})
    return (a[1] in b) || (a[2] in b)
end

function acc!(dst::Dict{PlaneKeyTuple2, Vector{PlaneKeyTuple2}}, a::PlaneKeyTuple2, b::PlaneKeyTuple2)
    if !(haskey(dst, a))
        dst[a] = Vector{PlaneKeyTuple2}()
    end
    push!(dst[a], b)
end

function compute_corner_loop(corners::Vector{PlaneKeyTuple2})::Vector{PlaneKeyTuple2}
    n = length(corners)
    if n < 3
        return []
    end
    
    neighbors = Dict{PlaneKeyTuple2, Vector{PlaneKeyTuple2}}()
    for a in corners
        for b in corners
            if a < b && are_connected(a, b)
                acc!(neighbors, a, b)
                acc!(neighbors, b, a)
            end
        end
    end

    visited = Set{PlaneKeyTuple2}()
    dst = Vector{PlaneKeyTuple2}()

    function visit!(x::PlaneKeyTuple2)
        @assert !(x in visited)
        push!(dst, x)
        push!(visited, x)
    end

    start = corners[1]
    at = start
    while true
        visit!(at)

        found = at
        for next in neighbors[at]
            if !(next in visited)
                found = next
                break
            end
        end

        if found == at
            break
        else
            at = found
        end
    end
    if length(visited) == n && are_connected(start, at)
        return dst
    else
        return []
    end
end

struct Vertex
    position::Vector{Float64}
    color::RgbColor
end

Facet = Tuple{Integer, Integer, Integer}

struct TriMesh
    vertices::Vector{Vertex}
    facets::Vector{Facet}
end

function mesh_from_physical_object(beam::Beam)
    polyhedron = beam.polyhedron
    vertex_count = length(polyhedron.corners)
    vertices = Vector{Vertex}()
    vertex_index_map = Dict{PlaneKeyTuple3, Integer}()

    for (k, v) in polyhedron.corners
        index = length(vertex_index_map)+1
        vertex_index_map[k] = index
        push!(vertices, Vertex(v, beam.specs.color))
    end

    facets = Vector{Facet}()
    for (k, v) in polyhedron.planes
        ks = plane_corner_keys(polyhedron, k)
        loop = compute_corner_loop(ks)

        function vind(i::Integer)
            (a, b) = loop[i]
            return vertex_index_map[ordered_triplet(k, a, b)]
        end
        
        m = length(loop)
        if 0 < m
            for j in 2:(m-1)
                facet = (vind(1), vind(j), vind(j+1))
                (i0, i1, i2) = facet
                (v0, v1, v2) = [vertices[i].position for i in facet]
                normal = cross(v1 - v0, v2 - v0)
                if dot(normal, v.normal) > 0
                    facet = reverse(facet)
                end
                push!(facets, facet)
            end
        end
    end
    return TriMesh(vertices, facets)
end

function wavefront_obj_string(mesh::TriMesh)
    buf = IOBuffer()
    for v in mesh.vertices
        (x, y, z) = v.position
        println(buf, string("v ", x, " ", y, " ", z))
    end
    for (i, j, k) in mesh.facets
        println(buf, string("f ", i, " ", j, " ", k))
    end
    return String(take!(buf))
end

struct CornerPosition
    key::PlaneKeyTuple3
    position::Vector{Float64}
end

# All coordinates are in a local coordinate system.
# The y axis points in the direction of the beam.
# The x axis points to the side.
# The z (which is not used) is in the normal of the plane.
struct BeamCuttingPlan
    # The plane that is used for the plan
    plane_key::PlaneKey

    # The sequence of corners
    corners::Vector{CornerPosition}

    # Any annotationss
    annotations::Vector{Annotation}

    # The bounding box
    bbox::BBox{Float64}
end

function plane_cog(polyhedron::Polyhedron, pk::PlaneKey)
    corners = [v for (k, v) in polyhedron.corners if pk in k]
    return (1.0/length(corners))*reduce(+, corners)
end

function plane_cog(beam::Beam, k::PlaneKey)
    return plane_cog(beam.polyhedron, k)
end

# Constructs a cutting plan for a beam.
function beam_cutting_plan(plane_key::PlaneKey, corners::Vector{CornerPosition}, annotations::Vector{Annotation})
    return BeamCuttingPlan(plane_key, corners, annotations, compute_bbox([corner.position for corner in corners]))
end

function cutting_plan(
    # The beam for which to generate a plane
    beam::Beam,

    # The plane for the plan
    k::PlaneKey,

    # The direction of the beam
    specified_beam_dir::Vector{Float64})

    # The normal of the plane is pointing inward.
    # Flip it so that it points outwards.
    z = -normalize(beam.polyhedron.planes[k].normal)
    
    x = cross(specified_beam_dir, z)
    @assert 0 < norm(x)
    x = normalize(x)

    y = cross(z, x)

    cog = plane_cog(beam, k)
    basis = [x y z]
    world_local = RigidTransform(basis, cog)
    local_world = invert(world_local)

    println(string("DET: ", det(basis)))

    loop = [ordered_triplet(pair, k) for pair in compute_corner_loop(plane_corner_keys(beam.polyhedron, k))]
    @assert 0 < length(loop)

    function localize(X)
        return transform_position(local_world, X)[1:2]
    end
    
    corners = [CornerPosition(corner, localize(beam.polyhedron.corners[corner])) for corner in loop]
    annotations = [@set annotation.position = localize(annotation.position) for annotation in get(beam.annotations, k, PersistentVector{Annotation}())]
    return beam_cutting_plan(k, corners, annotations)
end

function cutting_plan(beam::Beam, k::PlaneKey)
    return cutting_plan(beam, k, beam_dir(beam))
end

function plan_width(cp::BeamCuttingPlan)
    return width(cp.bbox.intervals[1])
end

function plan_length(cp::BeamCuttingPlan)
    return width(cp.bbox.intervals[2])
end

function lx_point(v)
    return lx.Point(v[1], v[2])
end

struct RenderConfig
    render_factor::Number
    fontsize::Number
end

const default_render_config = RenderConfig(100, 12)

function solve_plane_key(ikey::PlaneKeyTuple3, jkey::PlaneKeyTuple3)
    return none
end

function average(positions)
    return (1.0/length(positions))*reduce(+, positions)
end


function render_cutting_plan(cp::BeamCuttingPlan, render_config::RenderConfig)
    function local_pt(x)
        return lx_point(render_config.render_factor*x)
    end
    
    @lx.pdf begin
        lx.fontsize(render_config.fontsize)
        n = length(cp.corners)
        lx.setdash("solid")
        positions = [local_pt(corner.position) for corner in cp.corners]
        mid = average(positions)
        lx.poly(positions, :stroke, close=true)
        lx.label.([string(corner.key) for corner in cp.corners], lx.slope.(mid, positions), positions, offset=5)
        lx.rulers()
    end
end
# Rendering plans
# http://juliagraphics.github.io/Luxor.jl/stable/

function demo()
    bs = BeamSpecs(1.0, 3.0, default_beam_color)

    # Create a new beam that points in the X direction.
    beam = orient_beam(new_beam(bs), [1.0, 0.0, 0.0], local_y_dir)

    # Cut the beam at 0.0 and 4.5.
    a = NamedPlane(:a, plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    b = NamedPlane(:b, plane_at_pos([-1.0, 0.0, 0.0], [4.5, 0.0, 0.0]))

    beam = cut(a, cut(b, beam))

    
    # This is the plane key to render. Take the side of the beam. It should have length 4.5 and height 3.0
    k = :beam_X_lower

    dpspecs = DrillingPlaneSpecs(0.25, 2)
    drilling_dir = [0.0, 1.0, 0.0]
    beam_planes = generate_drilling_planes(beam, dpspecs, drilling_dir)

    println(beam_planes)

    

    # Make a plan for that plane
    plan = cutting_plan(beam, k)

    # Render the plan
    render_cutting_plan(plan, default_render_config)
end


export demo # Load the module and call Blueprint.demo() in the REPL.

end # module



# ONLY FOR DEBUG
#Blueprint.demo()
