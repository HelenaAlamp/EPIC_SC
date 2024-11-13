
struct Lattice <: AbstractVector{AbstractElement}
    # Lattice parameters
    sequence:: Vector{AbstractElement}
    Lattice() = new(Vector{AbstractElement}())
end


Base.size(lattice::Lattice) = size(lattice.sequence)
Base.length(lattice::Lattice) = length(lattice.sequence)
Base.getindex(lattice::Lattice, i::Int) = lattice.sequence[i]
Base.push!(lattice::Lattice, element::AbstractElement) = push!(lattice.sequence, element)
Base.push!(lattice::Lattice, elements::Vector{AbstractElement}) = push!(lattice.sequence, elements...)
Base.push!(lattice::Lattice, elements::Lattice) = push!(lattice.sequence, elements.sequence...)
Base.pop!(lattice::Lattice) = pop!(lattice.sequence)

function Lattice(seq::AbstractVector{AbstractElement})
    lattice = Lattice()
    for element in seq
        push!(lattice, element)
    end
    return lattice
end

function track!(beam::BunchedBeam, lattice::Lattice)
    for element in lattice
        track!(beam, element)
    end
    return nothing
end