# https://towardsdatascience.com/how-to-test-your-software-with-julia-4050379a9f3

# Run tests either by
#  1. C-c C-b from emacs
#  2. Command line: julia test/test_blueprint.jl

using Test
include("../src/blueprint.jl") # Alternative 1
#using blueprint               # Alternative 2
using FunctionalCollections
using LinearAlgebra

# Not needed:
#import .blueprint


@testset "Plane tests" begin
    @test blueprint.plane_at_pos([1, 2, 3], [3, 3, 7]).offset == 30.0

    plane = blueprint.plane_at_pos([1, 2, 3], [10, 11, 12])

    @test plane.normal[1] == 1.0
    @test plane.normal[2] == 2.0
    @test plane.offset == 68.0
    @test blueprint.evaluate(plane, [10, 11, 12]) == 0.0
    @test blueprint.evaluate(plane, [11, 13, 15]) == 14.0
    @test blueprint.evaluate(plane, [9, 9, 9]) == -14.0
end

@testset "Plane intersection" begin
    bp = blueprint
    a = bp.Plane([1.0, 0.0, 0.0], 0.0)
    b = bp.Plane([0.0, 1.0, 0.0], 0.0)
    line = bp.intersect(a, b)
    @test norm(line.dir - [0, 0, 1.0]) < 1.0e-6
    @test norm(line.pos - [0, 0, 0]) < 1.0e-6
end

@testset "Plane intersection 2" begin
    bp = blueprint
    a = bp.plane_at_pos([-1.0, 2.0, 0.0], [0.0, 0.0, 0.0])
    b = bp.plane_at_pos([10.0, 0.0, 0.0], [2.0, 0.0, 0.0])
    line = bp.intersect(b, a)
    @test line.dir[1] == 0.0
    @test line.dir[2] == 0.0
    @test 0 < line.dir[3]
    @test norm(line.pos - [2.0, 1.0, 0.0]) < 1.0e-6
end

@testset "Plane shadowing" begin
    bp = blueprint    
    @test bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 3.0]),
                         bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]))
    @test !bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.001], [0.0, 0.0, 3.0]),
                          bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]))
    @test !bp.shadowed_by(bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 4.0]),
                          bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 3.0]))
end

@testset "Polyhedron tests" begin
    bp = blueprint
    planes = @Persistent Dict(:a => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :b => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :c => bp.plane_at_pos([-1.0, -1.0, 0.0], [0.5, 0.5, 0.0]))
    polyhedron = bp.polyhedron_from_planes(planes)

    @test 3 == length(polyhedron.planes)
    @test 3 == length(polyhedron.bounded_lines)
    @test 0 == length(polyhedron.corners)
end

@testset "Polyhedron tests 2" begin
    bp = blueprint
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

@testset "Beam tests" begin
    specs = blueprint.BeamSpecs(1, 3)
    f = blueprint.beam_factory("Mjao", specs)
    beam = blueprint.new_beam!(f)
    @test beam.name == "Mjao0"
    @test f.counter == 1
end

@testset "Half-space test" begin
    bp = blueprint
    plane = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.5])
    @test bp.inside_halfspace(plane, [100.0, -2220.0, 0.6])
    @test !(bp.inside_halfspace(plane, [100.0, -2220.0, 0.4]))
end

@testset "Ordered pair test" begin
    @test blueprint.ordered_pair(:a, :b) == (:a, :b)
    @test blueprint.ordered_pair(:b, :a) == (:a, :b)
end

@testset "Test plane/line intersection" begin
    bp = blueprint
    @test bp.intersect(bp.plane_at_pos([0.0, 0.0, -1.0], [0.0, 0.0, 3.5]),
                       bp.ParameterizedLine([0.0, 0.0, 0.5], [0.0, 0.0, 0.0])).lambda == 7.0
    @test !bp.exists(
        bp.intersect(bp.plane_at_pos([0.0, 0.0, -1.0], [0.0, 0.0, 3.5]),
                     bp.ParameterizedLine([0.0, 3.0, 0.0], [0.0, 0.0, 0.0])))

end

@testset "Update line bounds test" begin
    bp = blueprint

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
