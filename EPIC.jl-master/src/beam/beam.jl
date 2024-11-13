
include("particle.jl")


struct BunchedBeam<: AbstractBeam 
    particle::ParticleType
    num_particle::Float64
    total_energy::Float64  
    momentum::Float64  # Design Momentum of the beam
    gamma::Float64  # Relativistic gamma
    beta::Float64  # Relativistic beta v/c
    num_macro::Int64 # Number of macro particles
    dist::StructArray{ps6d{Float64}} # 6D distribution
    inzindex::Vector{Int64} # in which z slice

    temp1::Vector{Float64} # temporary variable for track!
    temp3::Vector{Float64} # temporary variable for track!
    temp5::Vector{Float64} # temporary variable for track!
    temp2::Vector{Float64} # temporary variable for track!
    temp4::Vector{Float64} # temporary variable for track!


    emittance::Vector{Float64} # emittance in x, y, z
    centroid::Vector{Float64} # centroid in x, px, y, py, z, pz
    beamsize::Vector{Float64} # beamsize in x, px, y, py, z, pz
    moment2nd::Matrix{Float64} # 2nd moment matrix

    znbin::Int64 # number of z bins
    zhist::Vector{Float64} # z histogram for wakefield and spacecharge
    zhist_edges::Vector{Float64} # z histogram edges

    ztemp1::Vector{Float64} # temporary variable for track!
    ztemp2::Vector{Float64} # temporary variable for track!
    ztemp3::Vector{Float64} # temporary variable for track!
    ztemp4::Vector{Float64} # temporary variable for track!

    function BunchedBeam(particle::ParticleType, np::Float64, energy::Float64, nmacro::Int, 
                        emittance::Vector{Float64}, centroid::Vector{Float64}=zeros(6))
        momentum=sqrt(energy*energy-particle.mass*particle.mass)
        gamma=energy/particle.mass
        beta=momentum/energy
        znbin=Int64(round(sqrt(nmacro)))÷2
        new(particle, np, energy,momentum,gamma,beta,nmacro,
            StructArray{ps6d{Float64}}(undef,nmacro), Vector{Int64}(undef, nmacro), 
            Vector{Float64}(undef, nmacro), Vector{Float64}(undef, nmacro), Vector{Float64}(undef, nmacro), # temporary variable for track!
            Vector{Float64}(undef, nmacro), Vector{Float64}(undef, nmacro), # temporary variable for track!
            emittance, centroid, zeros(6), zeros(6,6),
            znbin, zeros(znbin), Vector{Float64}(undef, znbin+1), 
            Vector{Float64}(undef, znbin), Vector{Float64}(undef, znbin), Vector{Float64}(undef, znbin), Vector{Float64}(undef, znbin+1))
    end
end

#r0=e^2/(4*pi*eps_0*m_0*c^2)
include("dist.jl")


