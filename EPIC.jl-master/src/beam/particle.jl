struct  ParticleType  # Properties of charged particle
    charge::Float64
    mass::Float64
    atomnum::Float64
    classrad0::Float64
    radconst::Float64
    function ParticleType(charge, mass, atn)
        classrad0=charge*charge/(atn*mass)/4/π/55.26349406*1e-6
        radconst=4*π/3*classrad0/mass/mass/mass
        new(charge, mass, atn, classrad0, radconst)
    end
end

ParticleType(charge, mass)=ParticleType(charge, mass, 1.0)

const ELECTRON=ParticleType(-1.0, 0.51099895e6)
const PROTON=ParticleType(1.0, 938.27208816e6)
const GOLDION=ParticleType(79.0, 931.49410242e6, 197)

struct ps6d{T}  <: FieldVector{6, T} # 6D phase space
    x::T
    px::T
    y::T
    py::T
    z::T
    dp::T
end


