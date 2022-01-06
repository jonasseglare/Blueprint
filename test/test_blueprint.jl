# https://towardsdatascience.com/how-to-test-your-software-with-julia-4050379a9f3

# Run tests either by
#  1. C-c C-b from emacs
#  2. Command line: julia test/test_blueprint.jl

using Test
include("../src/blueprint.jl")

import .blueprint

@testset "Plane tests" begin
    @test blueprint.plane_at_pos([1, 2, 3], [3, 3, 7]).offset == 30.0

    plane = blueprint.plane_at_pos([1, 2, 3], [10, 11, 12])

    @test plane.normal[1] == 1.0
    @test plane.normal[2] == 2.0
    @test plane.offset == 68.0
end

@testset "Beam tests" begin
    specs = blueprint.BeamSpecs(1, 3)
    f = blueprint.beam_factory("Mjao", specs)
    beam = blueprint.new_beam!(f)
    @test beam.name == "Mjao0"
end
