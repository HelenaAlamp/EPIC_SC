#=
abstract type AbstractStrongBeamBeam <:AbstractElement end
using SpecialFunctions
using Statistics
using StaticArrays

struct StrongThinGaussianBeam <: AbstractStrongBeamBeam
    amplitude::Float64
    rmssizex::Float64
    rmssizey::Float64
    zloc::Float64
    xoffset::Float64
    yoffset::Float64
    StrongThinGaussianBeam(amp::Float64, rx::Float64, ry::Float64, zloc::Float64=0.0, xoff::Float64=0.0, yoff::Float64=0.0)=new(amp,rx,ry,zloc,xoff,yoff)
end

struct StrongGaussianBeam <: AbstractStrongBeamBeam  # Strong Beam with transverse Gaussian distribution
    particle::ParticleType
    num_particle::Float64  # Number of particles
    total_energy::Float64 # Total energy of the beam
    momentum::Float64  # Design Momentum of the beam
    gamma::Float64  # Relativistic gamma
    beta::Float64  # Relativistic beta v/c
    optics::AbstractOptics4D # optics @IP
    beamsize::Vector{Float64} # Beam size at IP
    nzslice::Int64 # Number of slices in z direction
    zslice_center::Vector{Float64} # z center of each slice
    zslice_npar::Vector{Float64} # amplitude of each slice
    function StrongGaussianBeam(particle::ParticleType, np::Float64, energy::Float64, op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)
        momentum=sqrt(energy*energy-particle.mass*particle.mass)  
        gamma=energy/particle.mass
        beta=momentum/energy
        new(particle,np,energy,momentum,gamma,beta, op, bs, Int64(nz), zeros(nz), zeros(nz))
    end  
end

function initilize_zslice!(beam::StrongGaussianBeam, profile::Symbol, slice_type::Symbol, zrange::Float64=5.0)
    zmin=-zrange*beam.beamsize[3]
    zmax=zrange*beam.beamsize[3]
    if profile == :gaussian
        if slice_type == :evenzsep
            zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
            beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
            beam.zslice_npar .= exp.(-0.5.*(beam.zslice_center.^2)./beam.beamsize[3]^2)
            beam.zslice_npar .= beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
        end
        if slice_type == :evennpar   # Here the zslicecenter is where split the npar in the slice, not the center of zposition
            npartedge=collect(range(0.0, stop=1.0, length=beam.nzslice+1))
            npartcenter=0.5.*(npartedge[1:end-1]+npartedge[2:end])  
            beam.zslice_center .= (sqrt(2.0)*beam.beamsize[3]).*erfinv.(2.0.*npartcenter.-1.0)
            beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
        end
    end

    if profile == :uniform
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
    end
    
end

function initilize_zslice!(beam::StrongGaussianBeam, zlist::Vector{Float64}, slice_type::Symbol) # Convert from distribution
    zmin=minimum(zlist)
    zmax=maximum(zlist)
    if slice_type == :evenzsep
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center=0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar=zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_npar[i]=sum((zlist.>=zedge[i]).&(zlist.<zedge[i+1]))
        end
        beam.zslice_npar=beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
    end
    if slice_type == :evennpar
        sort_zlist=sort(zlist)
        nzlist=length(sort_zlist)
        npartedge=Int.(floor.(collect(range(0, nzlist, length=beam.nzslice+1))))
        beam.zslice_npar .= (npartedge[2:end]-npartedge[1:end-1])/nzlist*beam.num_particle
        beam.zslice_center .= zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_center[i]=mean(sort_zlist[npartedge[i]+1:npartedge[i+1]])
        end
    end
end


function Bassetti_Erskine_xgty!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64) # x size greater than y
    # Only positive y is valid for this function
    # for y<0, Ex = Ex, Ey = -Ey
    if y < 0.0
        Bassetti_Erskine_xgty!(res, x, -y, σx, σy)
        res[2] = -res[2]
        nothing
        return
    end
    termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
	res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/σx/σy
    nothing
end

function Bassetti_Erskine_ygtx!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64) # x size greater than y
    # Only negative x is valid for this function
    # for x>0, Ex = -Ex, Ey = Ey
    if x > 0.0
        Bassetti_Erskine_ygtx!(res, -x, y, σx, σy)
        res[1] = -res[1]
        nothing
        return
    end
    termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
    res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/σx/σy
	nothing
end

function Bassetti_Erskine!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64)
    if σx > σy
        Bassetti_Erskine_xgty!(res, x, y, σx, σy)
        nothing
        return
    else
        Bassetti_Erskine_ygtx!(res, x, y, σx, σy)
        nothing
        return
    end
end

function track!(dist::AbstractVector{ps6d{T}}, temp1, temp2, temp3, temp4, temp5, sgb::StrongGaussianBeam, factor::Float64) where T
    #factor=wb.particle.classrad0/wb.gamma*wb.particle.charge*sgb.particle.charge = r_0/(gamma)*e*p
    
    lumi=0.0
    num_macro=length(dist.x)
    fieldvec = MVector{3,Float64}(0.0, 0.0, 0.0)

    @inbounds for i in 1:sgb.nzslice
        # temp1: collision zlocation, temp2: beamsize x, temp3: beamsize y, temp4: beta x, temp5: beta y
        temp1 .= (dist.z .+ sgb.zslice_center[i])./2.0
        temp4 .= sgb.optics.optics_x.beta .+ sgb.optics.optics_x.gamma .* temp1 .* temp1 .- 2.0 .* sgb.optics.optics_x.alpha .* temp1
        temp2 .= sgb.beamsize[1] .* sqrt.(temp4 ./ sgb.optics.optics_x.beta)
        temp5 .= sgb.optics.optics_y.beta .+ sgb.optics.optics_y.gamma .* temp1 .* temp1 .- 2.0 .* sgb.optics.optics_y.alpha .* temp1
        temp3 .= sgb.beamsize[2] .* sqrt.(temp5 ./ sgb.optics.optics_y.beta)
        
        # temp4 and temp5 are free to change now.

        dist.x .+= (dist.px .* temp1)
        dist.y .+= (dist.py .* temp1)
        dist.dp .-=  (dist.px .* dist.px .+ dist.py .* dist.py)./ 4.0
        
        slicelumi=0.0
        @inbounds for j in 1:num_macro
            Bassetti_Erskine!(fieldvec, dist.x[j], dist.y[j], temp2[j], temp3[j])
            dist.px[j] += (sgb.zslice_npar[i]*factor) * fieldvec[1]
            dist.py[j] += (sgb.zslice_npar[i]*factor) * fieldvec[2]
            slicelumi += fieldvec[3]
        end
       
        lumi += slicelumi * sgb.zslice_npar[i] #  Will do it outside* wb.num_particle / wb.num_macro
        
        dist.x .-= (dist.px .* temp1)
        dist.y .-= (dist.py .* temp1)
        dist.dp .+=  (dist.px .* dist.px .+ dist.py .* dist.py)./ 4.0

        
    end
    lumi

end

function track!(wb::BunchedBeam, sgb::StrongGaussianBeam)
    factor=wb.particle.classrad0/wb.gamma*wb.particle.charge*sgb.particle.charge
    lumi=track!(wb.dist, wb.temp1, wb.temp2, wb.temp3, wb.temp4, wb.temp5, sgb, factor)
    lumi *= wb.num_particle / wb.num_macro
end

=#