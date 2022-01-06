module blueprint

using LinearAlgebra

### Planes

## Represents the plane normal*X = offset
struct Plane{T}
    normal::Vector{T}
    offset::T
end

function plane_at_pos(normal::Vector{T}, pos::Vector{T}) where {T}
    Plane(normal, dot(pos, normal))
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

struct Beam
    name::String
    specs::BeamSpecs
    transform::RigidTransform{Float64}
end

function new_beam(name::String, specs::BeamSpecs)
    Beam(name, specs, identity_3d_transform)
end

function new_beam!(factory::BeamFactory)
    result = new_beam(string(factory.prefix, string(factory.counter)), factory.specs)
    factory.counter += 1
    result
end


# What should be included with `using`.
export plane_at_pos, BeamSpecs, beam_factory, new_beam

end # module
