# https://towardsdatascience.com/how-to-test-your-software-with-julia-4050379a9f3

# Run tests either by
#  1. C-c C-b from emacs
#  2. Command line: julia test/runtests.jl

using Test
using Setfield

include("../src/Blueprint.jl") # Alternative 1
#using Blueprint               # Alternative 2: Does not reload changes in REPL

using FunctionalCollections
using LinearAlgebra

# Not needed:
#import .Blueprint


@testset "Interval tests" begin
    bp = Blueprint
    x = bp.extend(bp.extend(bp.undefined_interval(), 3.0), 4.5)
    @test bp.is_defined(x)
    @test x.lower == 3.0
    @test x.upper == 4.5
end

@testset "Bounding box test" begin
    points = [[1.0, 2.3], [0.3, 5.9], [2.0, 0.1]]
    bp = Blueprint
    bbox = bp.compute_bbox(points)
    @test bbox.intervals[1].lower == 0.3
    @test bbox.intervals[1].upper == 2.0
    
    @test bbox.intervals[2].lower == 0.1
    @test bbox.intervals[2].upper == 5.9
end

@testset "Plane tests" begin
    @test Blueprint.plane_at_pos([1, 2, 3], [3, 3, 7]).offset == 30.0

    plane = Blueprint.plane_at_pos([1, 2, 3], [10, 11, 12])

    @test plane.normal[1] == 1.0
    @test plane.normal[2] == 2.0
    @test plane.offset == 68.0
    @test Blueprint.evaluate(plane, [10, 11, 12]) == 0.0
    @test Blueprint.evaluate(plane, [11, 13, 15]) == 14.0
    @test Blueprint.evaluate(plane, [9, 9, 9]) == -14.0
end

@testset "Plane intersection" begin
    bp = Blueprint
    a = bp.Plane([1.0, 0.0, 0.0], 0.0)
    b = bp.Plane([0.0, 1.0, 0.0], 0.0)
    line = bp.intersect(a, b)
    @test norm(line.dir - [0, 0, 1.0]) < 1.0e-6
    @test norm(line.pos - [0, 0, 0]) < 1.0e-6
end

@testset "Plane intersection 2" begin
    bp = Blueprint
    a = bp.plane_at_pos([-1.0, 2.0, 0.0], [0.0, 0.0, 0.0])
    b = bp.plane_at_pos([10.0, 0.0, 0.0], [2.0, 0.0, 0.0])
    line = bp.intersect(b, a)
    @test line.dir[1] == 0.0
    @test line.dir[2] == 0.0
    @test 0 < line.dir[3]
    @test norm(line.pos - [2.0, 1.0, 0.0]) < 1.0e-6
end

@testset "Plane shadowing" begin
    bp = Blueprint    
    @test bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 3.0]),
                         bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]))
    @test !bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.001], [0.0, 0.0, 3.0]),
                          bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]))
    @test !bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]),
                          bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 3.0]))
end

@testset "Parallel plane distance" begin
    bp = Blueprint

    a = bp.plane_at_pos([2.0, 0.0, 0.0], [1.0, 0.0, 0.0])
    b = bp.plane_at_pos([-100.0, 0.0, 0.0], [7.0, 0.0, 0.0])
    @test isapprox(6.0, bp.parallel_plane_distance(a, b), atol=1.0e-6)
    @test isapprox(6.0, bp.parallel_plane_distance(b, a), atol=1.0e-6)
end

@testset "Polyhedron tests" begin
    bp = Blueprint
    planes = @Persistent Dict(:a => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :b => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :c => bp.plane_at_pos([-1.0, -1.0, 0.0], [0.5, 0.5, 0.0]))
    polyhedron = bp.polyhedron_from_planes(planes)

    @test 3 == length(polyhedron.planes)
    @test 3 == length(polyhedron.bounded_lines)
    @test 0 == length(polyhedron.corners)
end

@testset "Polyhedron tests 2" begin
    bp = Blueprint
    planes = @Persistent Dict(:x => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :y => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :z => bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
                              :xyz => bp.plane_at_pos([-1.0, -1.0, -1.0], [1.0, 0.0, 0.0]))
                              
    polyhedron = bp.polyhedron_from_planes(planes)

    @test 4 == length(polyhedron.planes)
    @test 6 == length(polyhedron.bounded_lines)
    @test 4 == length(polyhedron.corners)

    @test [0.0, 0.0, 0.0] == polyhedron.corners[(:x, :y, :z)]
    @test [0.0, 0.0, 1.0] == polyhedron.corners[(:x, :xyz, :y)]
end

@testset "Polyhedron tests 3, redundant plane" begin
    bp = Blueprint
    planes = @Persistent Dict(:x => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :y => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :z => bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
                              :xyz => bp.plane_at_pos([-1.0, -1.0, -1.0], [1.0, 0.0, 0.0]),
                              :xyz2 => bp.plane_at_pos([-1.0, -1.0, -1.0], [1.1, 0.0, 0.0]))
                              
    polyhedron = bp.polyhedron_from_planes(planes)

    @test 4 == length(polyhedron.planes)
    @test 6 == length(polyhedron.bounded_lines)
    @test 4 == length(polyhedron.corners)

    @test [0.0, 0.0, 0.0] == polyhedron.corners[(:x, :y, :z)]
    @test [0.0, 0.0, 1.0] == polyhedron.corners[(:x, :xyz, :y)]
end

## Empty polyhedron

@testset "Empty polyhedron" begin
    bp = Blueprint
    planes = @Persistent Dict(:x => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :y => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :z => bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
                              :xyz => bp.plane_at_pos([-1.0, -1.0, -1.0], [-1.0, 0.0, 0.0]))
                              
    polyhedron = bp.polyhedron_from_planes(planes)

    @test !bp.has_corners(polyhedron)
    
    @test 4 == length(polyhedron.planes)
    @test 0 == length(polyhedron.bounded_lines)
    @test 0 == length(polyhedron.corners)
end


## Half space

@testset "Half-space test" begin
    bp = Blueprint
    plane = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.5])
    @test bp.inside_halfspace(plane, [100.0, -2220.0, 0.6])
    @test !(bp.inside_halfspace(plane, [100.0, -2220.0, 0.4]))
end

@testset "Ordered pair test" begin
    @test Blueprint.ordered_pair(:a, :b) == (:a, :b)
    @test Blueprint.ordered_pair(:b, :a) == (:a, :b)
end

@testset "Test plane/line intersection" begin
    bp = Blueprint
    @test bp.intersect(bp.plane_at_pos([0.0, 0.0, -1.0], [0.0, 0.0, 3.5]),
                       bp.ParameterizedLine([0.0, 0.0, 0.5], [0.0, 0.0, 0.0])).lambda == 7.0
    @test !bp.exists(
        bp.intersect(bp.plane_at_pos([0.0, 0.0, -1.0], [0.0, 0.0, 3.5]),
                     bp.ParameterizedLine([0.0, 3.0, 0.0], [0.0, 0.0, 0.0])))

end

@testset "Update line bounds test" begin
    bp = Blueprint

    ps = bp.default_polyhedron_settings()
    
    line = bp.ParameterizedLine([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
    bds = bp.initialize_line_bounds(line)

    A = 

    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([-1.0, 0.0, -1.0], [2.0, 0.0, 0.0]))

    @test bds.exists
    @test bds.lower == nothing
    @test bds.upper.value == 2.0

    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([-1.0, 0.0, -1.0], [3.0, 0.0, 0.0]))
    
    @test bds.exists
    @test bds.lower == nothing
    @test bds.upper.value == 2.0
    
    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([-1.0, 0.0, -1.0], [1.0, 0.0, 0.0]))
    
    @test bds.exists
    @test bds.lower == nothing
    @test bds.upper.value == 1.0

    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 1.0], [-4.5, 0.0, 0.0]))

    @test bds.exists
    @test bds.lower.value == -4.5
    @test bds.upper.value == 1.0

    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 1.0], [-4.6, 0.0, 0.0]))
    
    @test bds.exists
    @test bds.lower.value == -4.5
    @test bds.upper.value == 1.0
    
    bds = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 1.0], [-4.1, 0.0, 0.0]))
    
    @test bds.exists
    @test bds.lower.value == -4.1
    @test bds.upper.value == 1.0

    bds0 = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 1.0], [10.0, 0.0, 0.0]))

    @test !bds0.exists

    bds1 = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 0.0], [-10.0, 0.0, 0.0]))

    @test bds == bds1
    
    bds2 = bp.update_line_bounds(
        bds,
        :x,
        bp.plane_at_pos([1.0, 0.0, 0.0], [10.0, 0.0, 0.0]))

    @test !bds2.exists
end

@testset "Plane transforms" begin
    bp = Blueprint
    plane = bp.plane_at_pos([0.0, 1.0, 0.0], [2.0, 1.0, 0.0])

    # Try both rotation and translation, separately
    plane2 = bp.transform(bp.rigid_transform_from_xy_rotation(0.5*pi, 3), plane)
    @test isapprox(plane2.normal, [-1.0, 0.0, 0.0], atol=1.0e-6)
    @test isapprox(plane2.offset, 1.0, atol=1.0-6)

    plane3 = bp.transform(bp.rigid_transform_from_translation([0.5, 3.4, 0.0]), plane)
    @test plane3.normal == [0.0, 1.0, 0.0]
    @test plane3.offset == 4.4
end

@testset "Test add planes" begin
    bp = Blueprint
    base_poly = bp.polyhedron_from_planes(Dict(:x => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                                               :y => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                                               :z => bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])))
    @test 1 == length(base_poly.corners)
    polyhedron = bp.add_planes(base_poly, Dict(:xyz => bp.plane_at_pos([-1.0, -1.0, -1.0], [1.0, 0.0, 0.0])))
    @test 4 == length(polyhedron.corners)
    @test [0.0, 1.0, 0.0] == polyhedron.corners[(:x, :xyz, :z)]

    tpoly = bp.transform(bp.rigid_transform_from_xy_rotation(0.5*pi, 3), polyhedron)
    @test isapprox([-1.0, 0.0, 0.0], tpoly.corners[(:x, :xyz, :z)], atol=1.0e-6)
end

@testset "Test flip" begin
    bp = Blueprint
    plane = bp.plane_at_pos([3, 0.5, 0.25], [-3.0, 0.0, -4.0])
    plane2 = bp.flip(plane)

    X = [9.0, 4.7, 2.0]
    @test bp.evaluate(plane, X) == -bp.evaluate(plane2, X)
end

@testset "New beam test" begin
    bp = Blueprint
    specs = bp.beam_specs(1.0, 3.0)
    beam = bp.new_beam(specs)
    @test 4 == length(beam.polyhedron.planes)
    @test 4 == length(beam.polyhedron.bounded_lines)
    @test 0 == length(beam.polyhedron.corners)
end

@testset "Beam orientation transform" begin
    bp = Blueprint
    transform = bp.orient_beam_transform([1.0, 2.0, 0.0], [1.0, 1.0, 0.0])
    RtR = transform.rotation'*transform.rotation
    @test isapprox(RtR, Matrix(1.0I, 3, 3), atol=1.0e-6)

    z_world = bp.transform_direction(transform, [0.0, 0.0, 1.0])
    @test isapprox(z_world, normalize([1.0, 2.0, 0.0]), atol=1.0e-6)
end

@testset "Operations on transforms" begin
    bp = Blueprint
    a = bp.rigid_transform_from_xy_rotation(0.25pi, 3)
    b = bp.rigid_transform_from_translation([3.0, 2.0, 1.5])

    ab = bp.compose(a, b)

    X = [0.1, 0.2, 4.7]
    Y0 = bp.transform_position(a, bp.transform_position(b, X))
    Y1 = bp.transform_position(ab, X)

    ab_inv = bp.invert(ab)
    
    @test isapprox(Y0, Y1, atol=1.0e-6)
    @test isapprox(X, bp.transform_position(ab_inv, Y1), atol=1.0e-6)
end

@testset "Oriented beam test" begin
    bp = Blueprint
    specs = bp.beam_specs(1.0, 3.0)
    beam = bp.orient_beam(bp.new_beam(specs), [0.0, 1.0, 0.0], [1.0, 0.0, 0.0])
    @test !bp.has_corners(beam)
    @test isapprox(beam.polyhedron.planes[:beam_X_lower].normal, [0.0, 0.0, 1.0], atol=1.0e-6)
    @test isapprox(beam.polyhedron.planes[:beam_X_upper].normal, [0.0, 0.0, -1.0], atol=1.0e-6)
    @test isapprox(beam.polyhedron.planes[:beam_Y_lower].normal, [1.0, 0.0, 0.0], atol=1.0e-6)
    @test isapprox(beam.polyhedron.planes[:beam_Y_upper].normal, [-1.0, 0.0, 0.0], atol=1.0e-6)

    @test 4 == length(bp.bounding_points(beam))

    plane = bp.plane_at_pos([1.0, 0.0, 0.0], [-12.0, 0.0, 0.0])
    floor = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, -100.0])

    beam = bp.push_against(floor, bp.push_against(plane, beam))

    @test 0 == length(beam.polyhedron.corners)
    pts = bp.bounding_points(beam)

    @test 4 == length(pts)
    @test -12.0 == minimum(map(xyz -> xyz[1], pts))
    @test -100.0 == minimum(map(xyz -> xyz[3], pts))
    @test -9.0 == maximum(map(xyz -> xyz[1], pts))
    @test -99.0 == maximum(map(xyz -> xyz[3], pts))

    wall = bp.NamedPlane(:wall, bp.plane_at_pos([0.0, -1.0, 0.0], [0.0, 0.0, 0.0]))

    beam2 = bp.cut(wall, beam)
    @test bp.has_corners(beam2)
    @test 4 == length(beam2.polyhedron.corners)
end

@testset "Drill test" begin
    bp = Blueprint
    bs = bp.beam_specs(1.0, 3.0)

    A = bp.transform(bp.rigid_transform_from_translation([0.0, 0.0, 1.0]),
                     bp.orient_beam(bp.new_beam(bs), [1.0, 0.0, 0.0], bp.local_y_dir))
    B = bp.orient_beam(bp.new_beam(bs), [0.0, 1.0, 0.0], bp.local_y_dir)

    @test !bp.has_corners(A)
    
    Apt = bp.mid_point(A)
    Bpt = bp.mid_point(B)

    dir = bp.compute_drilling_direction(A, B)
    @test [0.0, 0.0, -1.0] == dir
    @test [0.0, 0.0, 1.0] == bp.compute_drilling_direction(
        A, bp.transform(bp.rigid_transform_from_translation([0.0, 0.0, 3.0]), B))
    @test [0.0, 0.0, -1.0] == bp.compute_drilling_direction(
        A, bp.transform(bp.rigid_transform_from_translation([0.0, 0.0, 0.3]), B))
    @test [0.0, 0.0, 1.0] == bp.compute_drilling_direction(B, A)
    
    dp_specs = bp.DrillingPlaneSpecs(0.25, 2)

    A_planes = bp.generate_drilling_planes(A, dp_specs, dir)
    B_planes = bp.generate_drilling_planes(B, dp_specs, dir)

    drills = bp.generate_drills(dir, A_planes, B_planes, bp.DrillSpecs(0))

    @test 4 == length(drills)
    
    for drill in drills
        @test drill.line.dir == dir
    end

    A = bp.drill(A, drills)
    @test 1 == length(A.annotations)
    drill_holes = A.annotations[bp.beam_Y_upper]
    @test 4 == length(drill_holes)
    for annotation in drill_holes
        (x, y, z) = annotation.position
        @test z == 4.0
    end

    A2 = bp.transform(bp.rigid_transform_from_translation([0.0, 1000.0, 0.0]), A)
    @test 1 == length(A2.annotations)
    @test 4 == length(A2.annotations[bp.beam_Y_upper])
    A2 = bp.drill(A2, drills)
    @test 1 == length(A2.annotations)
    @test 4 == length(A2.annotations[bp.beam_Y_upper])
    
    
    
    A3 = bp.transform(bp.rigid_transform_from_translation([0.0, 0.0001, 0.0]), A)
    @test 1 == length(A3.annotations)
    @test 4 == length(A3.annotations[bp.beam_Y_upper])
    A3 = bp.drill(A3, drills)
    @test 1 == length(A3.annotations)
    @test 8 == length(A3.annotations[bp.beam_Y_upper])
end

@testset "Obj export test" begin
    bp = Blueprint
    bs = bp.beam_specs(1.0, 3.0)
    beam = bp.orient_beam(bp.new_beam(bs), [1.0, 0.0, 0.0], bp.local_y_dir)

    a = bp.NamedPlane(:a, bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    b = bp.NamedPlane(:b, bp.plane_at_pos([-1.0, 0.0, 0.0], [4.5, 0.0, 0.0]))

    beam = bp.cut(a, bp.cut(b, beam))

    # Export wavefront obj
    mesh = bp.mesh_from_physical_object(beam)
    obj = bp.wavefront_obj_string(mesh)
    @test occursin("f 1 3 5", obj)

    
    # k = :beam_X_lower
    # @test haskey(beam.polyhedron.planes, k)

    # plan = bp.cutting_plan(beam, k)

    # @test bp.plan_width(plan) == 3.0
    # @test bp.plan_length(plan) == 4.5

    # @test 4 == length(plan.corners)

end

@testset "Cutting plan test" begin
    bp = Blueprint
    bs = bp.beam_specs(1.0, 3.0)

    # Create a new beam that points in the X direction.
    beam = bp.orient_beam(bp.new_beam(bs), [1.0, 0.0, 0.0], bp.local_y_dir)

    # Cut the beam at 0.0 and 4.5.
    a = bp.NamedPlane(:a, bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    b = bp.NamedPlane(:b, bp.plane_at_pos([-1.0, 0.0, 0.0], [4.5, 0.0, 0.0]))

    beam = bp.cut(a, bp.cut(b, beam))

    
    # This is the plane key to render. Take the side of the beam. It should have length 4.5 and height 3.0
    k = :beam_X_lower

    dpspecs = bp.DrillingPlaneSpecs(0.25, 2)
    drilling_dir = [0.0, 1.0, 0.0]
    
    beam_planes = bp.generate_drilling_planes(beam, dpspecs, drilling_dir)
    cut_planes = bp.generate_drilling_planes(bp.mid_point(beam), a.plane, b.plane, dpspecs)
    drills = bp.generate_drills(drilling_dir, beam_planes, cut_planes, bp.DrillSpecs(0))

    beam = bp.drill(beam, drills)
    

    # Make a plan for that plane
    plan = bp.cutting_plan(beam, k)

    @test 4 == length(plan.corners)
    @test 4 == length(plan.annotations)
end

@testset "Unique index test" begin
    bp = Blueprint
    m = Dict{bp.LabelSpec, Integer}()
    @test bp.generate_unique_index(m, bp.drill_label_spec) == 0
    @test bp.generate_unique_index(m, bp.drill_label_spec) == 1
    @test bp.generate_unique_index(m, bp.drill_label_spec) == 2
end

@testset "Project" begin
    bp = Blueprint
    plane = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.5])
    projected = bp.project(plane, [3.0, 4.5, 100.0])
    @test isapprox(projected, [3.0, 4.5, 0.5])
end

@testset "loop bounding planes test" begin
    bp = Blueprint
    base_loop = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

    for loop in [base_loop, reverse(base_loop)]
        planes = bp.loop_bounding_planes(loop)
        sample_points = [([-0.1, -0.1], false), ([0.1, 0.1], true)]
        function contained(x)
            for p in planes
                if bp.evaluate(p, x) < 0
                    return false
                end
            end
            return true
        end

        for (pos, expected) in sample_points
            @test expected == contained(pos)
        end
    end
end

@testset "Push against test" begin
    bp = Blueprint
    a = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    b = [[2.25, -1], [3.4, -1], [3.5, 10], [2.25, 10]]
    c = [[2.25, 5], [3.4, 5], [3.5, 10], [2.25, 10]]
    @test 1.25 == bp.push_loop_against_loop(a, b, [-1.0, 0.0])
    @test 1.25 == bp.push_loop_against_loop(b, a, [1.0, 0.0])
    @test 2.5 == bp.push_loop_against_loop(b, a, [0.5, 0.0])
    @test nothing == bp.push_loop_against_loop(a, c, [-1.0, 0.0])
end

## Simplification

@testset "Simplify edge multiples" begin
    bp = Blueprint
    src = [bp.IntervalEdge(0.3, 1),
           bp.IntervalEdge(0.6, -1),
           bp.IntervalEdge(0.6, 1),
           bp.IntervalEdge(0.75, -1)]
    dst = bp.simplify_edge_multiples(src)
    @test 2 == length(dst)
    @test dst[1].position == 0.3
    @test dst[2].position == 0.75
end

## Optimize it

@testset "Cut plan optimization" begin
    bp = Blueprint
    plans = [bp.sample_bcp(1.0, 0), bp.sample_bcp(2.0, 1), bp.sample_bcp(3.0, 2)]

    out = bp.pack(plans, 9, 0.25)

    @test 2 == length(out)
    @test 2 == length(out[1].plans)
    @test 1 == length(out[2].plans)
    @test 5 == bp.plan_length(out[1].plans[1])
    @test 4 == bp.plan_length(out[1].plans[2])
    @test 5.25 == bp.right(out[1].plans[1])
    @test 8.5 == bp.right(out[1].plans[2])
    @test 9.0 == bp.width(bp.bbox(out[1]).intervals[1])

    layout = out[1]
    edges = bp.generate_interval_edges(layout)

    intervals = bp.generate_incompressible_intervals(bp.simplify_edge_multiples(edges))
    @test 3 == length(intervals)

    xy = bp.generate_compression_function_points(intervals, 1.0)
    @test 6 == length(xy)
    @test (8.5, 5.5) == last(xy)

    @test 0 == length(bp.generate_incompressible_intervals(Vector{bp.IntervalEdge}()))
    @test 0 == length(bp.generate_compression_function_points(
        Vector{bp.DefinedInterval{Float64}}(), 1.0))

    @test 0.2 == bp.evaluate_piecewise_linear(xy, 0.2)
    @test 6.0 == bp.evaluate_piecewise_linear(xy, 9.0)
    @test 1.75 == bp.evaluate_piecewise_linear(xy, 2.75)

    f = bp.compression_function(layout, 1.0)
    @test 1.75 == f(2.75)
end

## Beam

@testset "Beam grouping and operations" begin
    bp = Blueprint
    bs = @set bp.beam_specs(1.0, 3.0).length = 6.0
    
    beam = bp.orient_beam(bp.new_beam(bs), [1.0, 0.0, 0.0], bp.local_y_dir)

    a = bp.NamedPlane(:a, bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
    b = bp.NamedPlane(:b, bp.plane_at_pos([-1.0, 0.0, 0.0], [4.5, 0.0, 0.0]))

    beam0 = bp.cut(a, bp.cut(b, beam))
    beam1 = bp.transform(bp.rigid_transform_from_translation([5.0, 0.0, 0.0]), beam0)

    @test 1 == length(bp.group([beam0]).components)
    @test 2 == length(bp.group([beam0, [[[[beam1]]]]]).components)
    @test 2 == length(bp.group(bp.group([beam0, [[[[beam1]]]]])).components)


    beams = bp.group([beam0, beam1])
    @test 2 == length(beams.components)
    @test 0.0 == bp.min_projection(a.plane, beams)
    @test -5.0 == bp.min_projection(b.plane, beams)

    beams2 = bp.transform(bp.rigid_transform_from_translation([2.5, 0.0, 0.0]), beams)
    @test 2.5 == bp.min_projection(a.plane, beams2)
    @test -7.5 == bp.min_projection(b.plane, beams2)

    beams3 = bp.push_against(b.plane, beams2)
    @test -5.0 == bp.min_projection(a.plane, beams3)
    @test 0.0 == bp.min_projection(b.plane, beams3)

    beamvec = bp.get_beams(beams3)

    @test beamvec[1].index == 0
    @test beamvec[2].index == 1

    groups = bp.group_by_specs(beamvec)
    @test 1 == length(groups)

    k = collect(keys(groups))[1]
    pgroup = groups[k]
    @test 2 == length(pgroup)

    packed_layouts = bp.pack_plans(groups)
    @test 1 == length(packed_layouts)

    layouts = packed_layouts[k]
    @test 2 == length(layouts)
end

## Numeric

@testset "Numeric string in base" begin
    base = collect("abc")
    bp = Blueprint
    @test "a" == bp.numeric_string_in_base(base, 0)
    @test "b" == bp.numeric_string_in_base(base, 1)
    @test "c" == bp.numeric_string_in_base(base, 2)
    @test "ba" == bp.numeric_string_in_base(base, 3)
end
