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
const Vector64 = Vector{Float64}

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

function bounding_points(bbox::BBox{T}) where {T}
    return [[ivl.lower for ivl in bbox.intervals], [ivl.upper for ivl in bbox.intervals]]
end

function compute_bbox(boxes::Vector{BBox{T}}) where {T}
    return compute_bbox([pt for bbox in boxes for pt in bounding_points(bbox)])
end


### Represents the plane normal*X = offset
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

function translate_normalized(plane::Plane{T}, amount::T) where {T}
    return translate(normalize_plane(plane), amount)
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


### Line

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

function rigid_transform_from_R(R::Matrix{T}) where {T}
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
        return dst
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

### Example:
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

    # Whether a side looks the same as the facing side
    flip_symmetric::Bool

    
    # The Z direction is the direction of the beam
end

default_beam_color = RgbColor(0.0, 0.0, 1.0)

function beam_specs(Xsize::Real, Ysize::Real)
    return BeamSpecs(Xsize, Ysize, default_beam_color, true)
end

function quadratic_beam_specs(size::Real)
    return beam_specs(size, size)
end

struct LabelSpec
    text::String
    short::String
end

struct Label
    spec::LabelSpec
    counter::Integer
end

function short_text(label::Label)
    return string(label.spec.short, label.counter)
end

abstract type AnnotationData end

struct Annotation
    position::Vector{Float64}
    label::Label
    data::AnnotationData
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
    T = rigid_transform_from_R(B*A')
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


function plane_interpolator(cog::Vector{Float64},
    lower_plane::Plane{Float64},
    upper_plane::Plane{Float64})

    lower_pos = project(lower_plane, cog)
    upper_pos = project(upper_plane, cog)

    lower_normal = lower_plane.normal
    upper_normal = sign(dot(upper_plane.normal, lower_plane.normal))*-1*upper_plane.normal

    function f(lambda)::Plane{Float64}
        lw = 1.0 - lambda
        uw = lambda
        pos = lw*lower_pos + uw*upper_pos
        normal = normalize(lw*lower_normal + uw*upper_normal)
        return plane_at_pos(normal, pos)
    end
    return f
end

function generate_drilling_planes(cog0::Vector{Float64},
    lower_plane0::Plane{Float64},
    upper_plane0::Plane{Float64},
    specs::DrillingPlaneSpecs)::Vector{Plane{Float64}}

    cog = 0.5*(project(lower_plane0, cog0) + project(upper_plane0, cog0))
    lower_plane = translate_normalized(lower_plane0, specs.marg)
    upper_plane = translate_normalized(upper_plane0, specs.marg)
    
    @assert dot(lower_plane.normal, upper_plane.normal) < 0
    if evaluate(lower_plane, cog) <= 0
        return Vector{Plane{Float64}}()
    end
    last = specs.count - 1
    step = 1.0/last
    f = plane_interpolator(cog, lower_plane, upper_plane)
    return map(i -> f(step*i), 0:last)
end
    


function generate_drilling_planes(
    beam::Beam, specs::DrillingPlaneSpecs, drilling_dir::Vector{Float64})::Vector{Plane{Float64}}
    dim = base_drilling_dim(beam, drilling_dir)
    (lower_key, upper_key) = plane_keys_per_beam_dim[dim]
    lower_plane = beam.polyhedron.planes[lower_key]
    upper_plane = beam.polyhedron.planes[upper_key]
    return generate_drilling_planes(mid_point(beam), lower_plane, upper_plane, specs)
end


function generate_unique_index(m::Dict{K, Integer}, k::K) where {K}
    i = get(m, k, 0)
    m[k] = i+1
    return i   
end

const drill_label_spec = LabelSpec("Drill", "D")

struct DrillSpecs <: AnnotationData
    radius::Number
end

struct Drilling
    label_spec::LabelSpec
    line::ParameterizedLine{Float64}
    specs::DrillSpecs
end

function generate_drills(drilling_dir::Vector{Float64},
    a_planes::Vector{Plane{Float64}},
    b_planes::Vector{Plane{Float64}},
    specs::DrillSpecs)::Vector{Drilling}
    
    dst = Vector{Drilling}()
    for a in a_planes
        for b in b_planes
            push!(dst, Drilling(drill_label_spec, direct_like(intersect(a, b), drilling_dir), specs))
        end
    end
    return dst
end

function drill(beam::Beam, drills::AbstractVector{Drilling})
    new_annotations = Dict{PlaneKey, PersistentVector{Annotation}}()
    for (plane_key, plane) in beam.polyhedron.planes
        plane = normalize_plane(plane)
        dst = get(beam.annotations, plane_key, pvec(Vector{Annotation}()))
        counters = Dict{LabelSpec, Integer}()
        for drill in drills
            if dot(drill.line.dir, plane.normal) > 0
                intersection = intersect(plane, drill.line)
                if exists(intersection)
                    counter = generate_unique_index(counters, drill.label_spec)
                    dst = push(dst, Annotation(evaluate(drill.line, intersection.lambda), Label(drill.label_spec, counter), drill.specs))
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

abstract type AbstractCuttingPlan end

# All coordinates are in a local coordinate system.
#
# When looking at the beam side from outside the beam:
#
# * The x axis points in the direction of the beam (right)
# * The y axis points to the side (down)
# * The z (which is not used) points away, in toward the inside of the polyhedron
#
#
struct BeamCuttingPlan <: AbstractCuttingPlan
    # The plane that is used for the plan
    plane_key::PlaneKey

    flip_symmetric::Bool
    
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
function beam_cutting_plan(plane_key::PlaneKey, flipsym::Bool, corners::Vector{CornerPosition}, annotations::Vector{Annotation})
    return BeamCuttingPlan(plane_key, flipsym, corners, annotations, compute_bbox([corner.position for corner in corners]))
end

function transform(rt::RigidTransform{Float64}, plan::BeamCuttingPlan)
    return beam_cutting_plan(
    plan.plane_key,
    plan.flip_symmetric,
    [@set corner.position = transform_position(rt, corner.position) for corner in plan.corners],
    [@set annotation.position = transform_position(rt, annotation.position) for annotation in plan.annotations])
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
    z = normalize(beam.polyhedron.planes[k].normal)
    
    y = cross(z, specified_beam_dir)
    @assert 0 < norm(y)
    y = normalize(y)

    x = cross(y, z)

    cog = plane_cog(beam, k)
    basis = [x y z]

    @assert abs(det(basis) - 1.0) < 1.0e-3
    
    world_local = RigidTransform(basis, cog)
    local_world = invert(world_local)

    

    loop = [ordered_triplet(pair, k) for pair in compute_corner_loop(plane_corner_keys(beam.polyhedron, k))]
    @assert 0 < length(loop)

    function localize(X)
        return transform_position(local_world, X)[1:2]
    end
    
    corners = [CornerPosition(corner, localize(beam.polyhedron.corners[corner])) for corner in loop]
    annotations = [@set annotation.position = localize(annotation.position) for annotation in get(beam.annotations, k, PersistentVector{Annotation}())]
    return beam_cutting_plan(k, beam.specs.flip_symmetric, corners, annotations)
end

function cutting_plan(beam::Beam, k::PlaneKey)
    return cutting_plan(beam, k, beam_dir(beam))
end

function lx_point(v)
    return lx.Point(v[1], v[2])
end

struct RenderConfig
    render_factor::Number
    fontsize::Number
    offset::Number
    marker_size::Number

    width::Number
    height::Number
    filename::String
end

const default_render_config = RenderConfig(100, 12, 5, 5, 600, 600, "/tmp/blueprint_sketch.pdf")

function solve_plane_key(ikey::PlaneKeyTuple3, jkey::PlaneKeyTuple3)
    return none
end

function average(positions)
    return (1.0/length(positions))*reduce(+, positions)
end


function render_pdf(body_fn, render_config::RenderConfig)
    @lx.pdf begin
        function local_pt(x)
            return lx_point(render_config.render_factor*x)
        end
        body_fn(local_pt)
    end render_config.width render_config.height render_config.filename
end

function render(cp::BeamCuttingPlan, render_config::RenderConfig)
    function sub(local_pt)
        @info "The matrix is" lx.getmatrix()
        offset = render_config.offset
        
        lx.fontsize(render_config.fontsize)
        n = length(cp.corners)
        lx.setdash("solid")
        positions = [local_pt(corner.position) for corner in cp.corners]
        mid = average(positions)
        lx.poly(positions, :stroke, close=true)
        lx.label.([string(corner.key) for corner in cp.corners], lx.slope.(mid, positions), positions, offset=offset)

        apos = [local_pt(annotation.position) for annotation in cp.annotations]

        lx.circle.(apos, render_config.marker_size, :fill)
        lx.label.([short_text(annotation.label) for annotation in cp.annotations], :SE, apos, offset=offset)

        lx.rulers()
     end
    render_pdf(sub, render_config)
end
# Rendering plans
# http://juliagraphics.github.io/Luxor.jl/stable/

function loop_bounding_planes(positions::Vector{Vector{Float64}})
    c = average(positions)
    n = length(positions)
    result = Vector{Plane{Float64}}()
    for i in 1:n
        src = positions[i]
        dst = positions[mod1(i + 1, n)]
        diff = dst - src
        normal = [-diff[2], diff[1]]
        plane = plane_at_pos(normal, src)
        if evaluate(plane, c) < 0
            plane = flip(plane)
        end
        push!(result, plane)
    end
    return result
end


function loop_bounding_planes(loop_corners::Vector{CornerPosition})
    return loop_bounding_planes(map(cp -> cp.position, loop_corners))
end

function solve_push_segments(left, pos, dir)
    (x, y) = left
    A = [x-y -dir]
    B = pos-y
    try
        (lambda, alpha) = A\B
        if 0 <= lambda && lambda <= 1.0
            return alpha
        else
            return nothing
        end
    catch e
        return nothing
    end
end

function push_loop_against_loop(
    corners_a::Vector{Vector64},
    corners_b::Vector{Vector64},
    direction::Vector64)

    planes_a = loop_bounding_planes(corners_a)
    planes_b = loop_bounding_planes(corners_b)

    an = length(corners_a)
    bn = length(corners_b)

    
    best::Union{Float64, Nothing} = nothing

    function get_segment(corners, index)
        return [corners[index], corners[mod1(index + 1, length(corners))]]
    end
    
    for i in 1:an
        adot = dot(direction, planes_a[i].normal)
        for j in 1:bn
            bdot = dot(direction, planes_b[j].normal)

            if adot*bdot <= 0

                src_dst = [get_segment(corners_a, i), get_segment(corners_b, j)]
                
                for k in 0:1
                    left = src_dst[1 + k]
                    dir = (1 - 2*k)*direction
                    (r0, r1) = src_dst[1 + (1 - k)]

                    x = dot(dir, r0) > dot(dir, r1) ? r0 : r1
                    
                    sol = solve_push_segments(left, x, dir)
                    if sol != nothing && (best == nothing || sol < best)
                        best = sol
                    end
                end
            end
        end
    end
    return best
end

function cutting_plan_loop(plan::BeamCuttingPlan)
    return [corner.position for corner in plan.corners]
end

function align_plan(plan::BeamCuttingPlan)
    offset = [plan.bbox.intervals[1].lower, plan.bbox.intervals[2].lower]
    return transform(rigid_transform_from_translation(-offset), plan)
end

function plan_length(plan::BeamCuttingPlan)
    return width(plan.bbox.intervals[1])
end

function longest_plan(plans::Vector{BeamCuttingPlan})
    return reduce((a, b) -> plan_length(a) > plan_length(b) ? a : b, plans)
end

function generate_transformed_plans(plan::BeamCuttingPlan)
    rotated = transform(rigid_transform_from_R([-1.0 0.0; 0.0 -1.0]), plan)
    M = rigid_transform_from_R([1.0 0.0; 0.0 -1.0])
    dst = [align_plan(plan), align_plan(rotated)]
    if plan.flip_symmetric
        push!(dst, align_plan(transform(M, plan)))
        push!(dst, align_plan(transform(M, rotated)))
    end
    return dst
    
end

function right(plan::BeamCuttingPlan)
    return plan.bbox.intervals[1].upper
end

struct BeamLayout
    beam_length::Number
    plans::Vector{BeamCuttingPlan}
end

function pack(plans::Vector{BeamCuttingPlan}, beam_length::Number, margin::Number)
    n = length(plans)
    plan_map = Dict{Int32, BeamCuttingPlan}()
    for plan in plans
        plan_map[length(plan_map)] = plan
    end

    result = Vector{BeamLayout}()    
    stop_loop = [[0.0, 0.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0]]
    dir = [-1.0, 0.0]
    
    current_beam = Vector{BeamCuttingPlan}()
    last_loop = stop_loop


    function push_layout()
        push!(result, BeamLayout(beam_length, current_beam))
        current_beam = Vector{BeamCuttingPlan}()
        last_loop = stop_loop
    end

    while 0 < length(plan_map)
        best_key = nothing
        best_plan = nothing
        best_measures = (0, -beam_length)
        for (k, plan) in plan_map
            genplans = generate_transformed_plans(plan)
            for gplan in genplans
                loop = cutting_plan_loop(gplan)
                amount = push_loop_against_loop(last_loop, loop, dir)
                adjusted_plan = transform(rigid_transform_from_translation((amount - margin)*dir), gplan)
                new_length = plan_length(adjusted_plan)
                new_right = right(adjusted_plan)
                new_measures = (new_length, -new_right)
                if new_right <= beam_length - margin && best_measures < new_measures
                    best_key = k
                    best_plan = adjusted_plan
                    best_measures = new_measures
                end
            end
        end

        if best_plan == nothing
            if length(current_beam) == 0
                error("Cutting plan does not fit in beam")
            else
                push_layout()
            end
        else
            push!(current_beam, best_plan)
            last_loop = cutting_plan_loop(best_plan)
            delete!(plan_map, best_key)
        end
    end

    if 0 < length(current_beam)
        push_layout()
    end
    
    return result
end

function bbox(beam_layout::BeamLayout)
    dst = compute_bbox([plan.bbox for plan in beam_layout.plans])
    @assert 0 <= dst.intervals[1].lower
    @assert dst.intervals[1].upper <= beam_layout.beam_length
    return BBox{Float64}([DefinedInterval{Float64}(0, beam_layout.beam_length), dst.intervals[2]])
end

struct IntervalEdge
    position::Float64
    step::Integer
end

function generate_interval_edges(plan::BeamCuttingPlan, tol=1.0e-6)
    positions = [c.position for c in plan.corners]
    n = length(positions)
    result = Vector{IntervalEdge}()
    for i in 1:n
        src = positions[i]
        dst = positions[mod1(i+1, n)]
        diff = dst - src
        dist = norm(diff)
        flat = dist > tol && abs(diff[2])/dist < tol
        if !flat
            (rising, falling) = sort([src[1], dst[1]])
            push!(result, IntervalEdge(rising, 1))
            push!(result, IntervalEdge(falling, -1))
        end
    end
    return result
end

function generate_interval_edges(layout::BeamLayout, tol=1.0e-6)
    edges = [edge for plan in layout.plans for edge in generate_interval_edges(plan, tol)]
    return sort(edges, by=x->x.position)
end

function simplify_edge_multiples(edges::Vector{IntervalEdge})
    dst0 = Vector{IntervalEdge}()
    if length(edges) == 0
        return dst0
    end
    
    last_pos = edges[1].position
    last_step = 0
    for edge in edges
        if edge.position != last_pos
            if last_step != 0
                push!(dst0, IntervalEdge(last_pos, last_step))
            end
            last_pos = edge.position
            last_step = edge.step
        else
            last_step += edge.step
        end
    end
    if last_step != 0
        push!(dst0, IntervalEdge(last_pos, last_step))
    end
    return dst0
end

function generate_incompressible_intervals(edges_no_multiples::Vector{IntervalEdge})
    current_level = 0
    start_at = 0
    dst = Vector{DefinedInterval{Float64}}()
    for edge in edges_no_multiples
        next_level = current_level + edge.step
        if current_level == 0 && next_level == 1
            start_at = edge.position
        elseif current_level == 1 && next_level == 0
            push!(dst, DefinedInterval{Float64}(start_at, edge.position))
        end
        current_level = next_level
    end
    return dst
end

function generate_compression_function_points(
    intervals::Vector{DefinedInterval{Float64}},
    max_step::Float64)
    
    dst = Vector{Tuple{Float64, Float64}}()
    n = length(intervals)
    if n == 0
        return dst
    end
    push!(dst, (intervals[1].lower, intervals[1].lower))
    push!(dst, (intervals[1].upper, intervals[1].upper))

    for i in 2:n
        left = intervals[i-1]
        right = intervals[i]
        (x, y0) = last(dst)
        step = right.lower - left.upper
        y1 = y0 + min(max_step, step)
        push!(dst, (right.lower, y1))
        push!(dst, (right.upper, y1 + (right.upper - right.lower)))
    end
    return dst
end

function evaluate_piecewise_linear_inner(xy_pairs::Vector{Tuple{Float64, Float64}},
                                         x::Float64, l::Int64, r::Int64)
    n = r - l
    @assert 1 <= n
    if n == 1
        (x0, y0) = xy_pairs[l]
        (x1, y1) = xy_pairs[r]
        @assert x0 < x1
        lambda = (x - x0)/(x1 - x0)
        return (1.0 - lambda)*y0 + lambda*y1
    else
        mid = l + div(n, 2)
        if x < xy_pairs[mid][1]
            return evaluate_piecewise_linear_inner(xy_pairs, x, l, mid)
        else
            return evaluate_piecewise_linear_inner(xy_pairs, x, mid, r)
        end
    end        
end
    

function evaluate_piecewise_linear(xy_pairs::Vector{Tuple{Float64, Float64}}, x)
    if length(xy_pairs) == 0
        return x
    end

    if x <= xy_pairs[1][1]
        return x
    end

    lxy = last(xy_pairs)
    if lxy[1] <= x
        return lxy[2] + (x - lxy[1])
    end

    return evaluate_piecewise_linear_inner(xy_pairs, x, 1, length(xy_pairs))
end

function compression_function(layout::BeamLayout, max_step::Float64, tol=1.0e-6)
    xy = generate_compression_function_points(
    generate_incompressible_intervals(
        simplify_edge_multiples(
            generate_interval_edges(layout, tol))), max_step)
    return x -> evaluate_piecewise_linear(xy, x)
end


######################## Samples code


function sample_bcp(len::Float64)
    bp = Blueprint
    k = :a
    corners = [bp.CornerPosition((:a, :b, :c), [2.0 + len, 0.0]),
               bp.CornerPosition((:a, :c, :d), [0.0, 0.0]),
               bp.CornerPosition((:a, :d, :e), [1.0, 1.0]),
               bp.CornerPosition((:a, :b, :e), [1.0 + len, 1.0])]
    annotations = Vector{bp.Annotation}()
    return bp.beam_cutting_plan(k, true, corners, annotations)
end

function demo0()
    bs = BeamSpecs(1.0, 2.0, default_beam_color, true)

    # Create a new beam that points in the X direction.
    beam = orient_beam(new_beam(bs), [1.0, 0.0, 0.0], local_y_dir)

    # Cut the beam at 0.0 and 4.5.
    a = NamedPlane(:a, plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    b = NamedPlane(:b, plane_at_pos([-1.0, 0.0, 0.3], [4.5, 0.0, 0.0]))

    beam = cut(a, cut(b, beam))

    
    # This is the plane key to render. Take the side of the beam. It should have length 4.5 and height 3.0
    k = :beam_X_lower

    dpspecs = DrillingPlaneSpecs(0.25, 2)
    drilling_dir = [0.0, 1.0, 0.0]
    
    beam_planes = generate_drilling_planes(beam, dpspecs, drilling_dir)
    cut_planes = translate_normalized.([a.plane, b.plane], dpspecs.marg)
    drills = generate_drills(drilling_dir, beam_planes, cut_planes, DrillSpecs(0.003))

    beam = drill(beam, drills)
    

    # Make a plan for that plane
    plan = cutting_plan(beam, k)

    println(string("Corners: ", length(plan.corners)))
    println(string("Annotations: ", length(plan.annotations)))

    # Render the plan
    render(plan, default_render_config)
end

function demo1()
    plans = [sample_bcp(1.0), sample_bcp(2.0), sample_bcp(3.0)]
    out = bp.pack(plans, 9, 0.25)

    render(out[1], default_render_config)
end


#export demo # Load the module and call Blueprint.demo() in the REPL.

end # module



# ONLY FOR DEBUG
#Blueprint.demo()
