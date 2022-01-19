# https://towardsdatascience.com/how-to-test-your-software-with-julia-4050379a9f3

# Run tests either by
#  1. C-c C-b from emacs
#  2. Command line: julia test/test_blueprint.jl

using Test
include("../src/blueprint.jl")

using FunctionalCollections
using LinearAlgebra

import .blueprint

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
    #if line == nothing
    #    print("Nothing")
    #end
end

@testset "Polyhedron tests" begin
    bp = blueprint
    planes = @Persistent Dict(:a => bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
                              :b => bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                              :c => bp.plane_at_pos([-1.0, -1.0, 0.0], [0.5, 0.5, 0.0]))
    polyhedron = bp.polyhedron_from_planes(planes)
    print(polyhedron)
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

@testset "Update line bounds test" begin
    bp = blueprint

    ps = bp.default_polyhedron_settings()
    
    line = bp.ParameterizedLine([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    bds = bp.initialize_line_bounds(line)

    A = bp.plane_at_pos([-1.0, 0.0, -1.0], [2.0, 0.0, 0.0])

    bds = bp.update_line_bounds(bds, A, ps.marg)

    println(string("BOUDNS: ", bds))
end
