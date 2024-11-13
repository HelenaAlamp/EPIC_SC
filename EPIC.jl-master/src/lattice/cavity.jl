

struct CrabCavity <: AbstractCavity
    f::Float64 # frequency
    v::Float64 # voltage
    k::Float64 # wavenumber
    ϕ::Float64 # phase
    errors::AbstractArray # 1: Voltage error, 2: Phase error
    CrabCavity(freq,nv, phase, error)=new(freq,nv,2*π*freq/2.99792458e8, phase, error)
end

CrabCavity(freq, nv, phase) = CrabCavity(freq, nv, phase, StaticArrays.@MVector [0.0,0.0])
CrabCavity(freq, nv) = CrabCavity(freq, nv, 0.0, StaticArrays.@MVector [0.0,0.0])




function track!(ps6dcoor::AbstractVector{ps6d{T}}, cc::CrabCavity, b2E, ang, sinang, cosang) where T
   ang .= (-cc.k) .* ps6dcoor.z .+ cc.ϕ
   sinang .= sin.(ang)
   cosang .= cos.(ang)
   ps6dcoor.px .+= (cc.v/b2E) .* sinang   #px   dx'=v*sin(-kz+phi0) = - dh/dx ==> h~ -v*sin(-kz+phi0)*x
   ps6dcoor.dp .+= (-cc.k* cc.v/b2E) .* ps6dcoor.x .* cosang  # delta dδ = -dh/dz = k v cos(kz) x
   return nothing
end

function track!(beam::AbstractBeam, cc::CrabCavity)
   β2E=beam.beta*beam.beta*beam.total_energy
   track!(beam.dist, cc, β2E, beam.temp1, beam.temp2, beam.temp3)
end

struct easyCrabCavity <: AbstractCavity
    f::Float64 # frequency
    halfθc::Float64 # voltage
    k::Float64 # wavenumber
    ϕ::Float64 # phase
    errors::AbstractArray # 1: Voltage error, 2: Phase error
    easyCrabCavity(freq, hθc, phase, error)=new(freq,hθc,2*π*freq/2.99792458e8, phase, error)
end
easyCrabCavity(freq, hθc, phase) = easyCrabCavity(freq, hθc, phase, StaticArrays.@MVector [0.0,0.0])
easyCrabCavity(freq, hθc) = easyCrabCavity(freq, hθc, 0.0, StaticArrays.@MVector [0.0,0.0])

function track!(ps6dcoor::AbstractVector{ps6d{T}}, cc::easyCrabCavity, ang, sinang, cosang) where T
   ang .= (-cc.k) .* ps6dcoor.z .+ cc.ϕ
   sinang .= sin.(ang)
   cosang .= cos.(ang)
   ps6dcoor.x .+= (cc.halfθc/cc.k) .* sinang    #px   dx=tc/k*sin(-kz+phi0) =  dh/dpx ==> h~ tc/k*sin(-kz+phi0)*px
   ps6dcoor.dp .+= (cc.halfθc) .* ps6dcoor.px .* cosang  # delta dδ = -dh/dz = tc cos(-kz+phi0) * px
   return nothing
end

function track!(beam::AbstractBeam, cc::easyCrabCavity)
   track!(beam.dist, cc, beam.temp1, beam.temp2, beam.temp3)
end




mutable struct AccelCavity <: AbstractCavity
   f::Float64 # frequency
   v::Float64 # voltage
   k::Float64 # wavenumber
   h::Int64 # harmonic number
   ϕs::Float64 # synchronous phase π/2 for accelerating on crest
   
   AccelCavity(freq,nv,h,ϕs)=new(freq,nv, 2*π*freq/2.99792458e8, h, ϕs)
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, ac::AccelCavity, β2E::Float64, sv) where T
   v_β2E=ac.v/β2E
   sv .= sin.((-ac.k) .* ps6dcoor.z .+ ac.ϕs) .- sin(ac.ϕs)
   ps6dcoor.dp .+= v_β2E .* sv   # update delta p
   return nothing
end



function track!(beam::AbstractBeam, ac::AccelCavity)
   β2E=beam.beta*beam.beta*beam.total_energy
   track!(beam.dist, ac, β2E, beam.temp1)
end



