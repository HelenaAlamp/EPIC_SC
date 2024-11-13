using LinearAlgebra
using Distributions

const distribution1D_type=Set([:custom, :gaussian, :uniform])
const zslice_type=Set([:evenzsep, :evennpar])

function get_centroid!(beam::BunchedBeam)
    # calculate centroid
    beam.centroid[1] = sum(beam.dist.x)/beam.num_macro
    beam.centroid[2] = sum(beam.dist.px)/beam.num_macro
    beam.centroid[3] = sum(beam.dist.y)/beam.num_macro
    beam.centroid[4] = sum(beam.dist.py)/beam.num_macro
    beam.centroid[5] = sum(beam.dist.z)/beam.num_macro
    beam.centroid[6] = sum(beam.dist.dp)/beam.num_macro
    return nothing
end


function get_2nd_moment!(beam::BunchedBeam)
    # calculate emittance
    
    ## Better way to do this?  Really messy but works...
    beam.moment2nd .= 0.0
    beam.temp1 .= beam.dist.x .* beam.dist.x;  beam.moment2nd[1,1] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.x .* beam.dist.px; beam.moment2nd[1,2] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.x .* beam.dist.y;  beam.moment2nd[1,3] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.x .* beam.dist.py; beam.moment2nd[1,4] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.x .* beam.dist.z;  beam.moment2nd[1,5] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.x .* beam.dist.dp; beam.moment2nd[1,6] = sum(beam.temp1)/beam.num_macro

    beam.temp1 .= beam.dist.px .* beam.dist.px; beam.moment2nd[2,2] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.px .* beam.dist.y;  beam.moment2nd[2,3] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.px .* beam.dist.py; beam.moment2nd[2,4] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.px .* beam.dist.z;  beam.moment2nd[2,5] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.px .* beam.dist.dp; beam.moment2nd[2,6] = sum(beam.temp1)/beam.num_macro

    beam.temp1 .= beam.dist.y .* beam.dist.y;  beam.moment2nd[3,3] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.y .* beam.dist.py; beam.moment2nd[3,4] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.y .* beam.dist.z;  beam.moment2nd[3,5] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.y .* beam.dist.dp; beam.moment2nd[3,6] = sum(beam.temp1)/beam.num_macro

    beam.temp1 .= beam.dist.py .* beam.dist.py; beam.moment2nd[4,4] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.py .* beam.dist.z;  beam.moment2nd[4,5] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.py .* beam.dist.dp; beam.moment2nd[4,6] = sum(beam.temp1)/beam.num_macro

    beam.temp1 .= beam.dist.z .* beam.dist.z;  beam.moment2nd[5,5] = sum(beam.temp1)/beam.num_macro
    beam.temp1 .= beam.dist.z .* beam.dist.dp; beam.moment2nd[5,6] = sum(beam.temp1)/beam.num_macro

    beam.temp1 .= beam.dist.dp .* beam.dist.dp; beam.moment2nd[6,6] = sum(beam.temp1)/beam.num_macro 

    
    beam.moment2nd .= beam.moment2nd+beam.moment2nd'-diagm(diag(beam.moment2nd))
    return nothing
end

function get_emittance!(beam::BunchedBeam)
    # calculate emittance
    get_centroid!(beam)
    get_2nd_moment!(beam)
    @inbounds for i in 1:6
        for j in 1:6
            beam.moment2nd[i,j] = beam.moment2nd[i,j] - beam.centroid[i]*beam.centroid[j]
        end
    end
    
    @inbounds for i in 1:3
        beam.emittance[i] = sqrt(det(@view beam.moment2nd[2*i-1:2*i,2*i-1:2*i]))
    end
    @inbounds for i in 1:6
        beam.beamsize[i] = sqrt(beam.moment2nd[i,i])
    end
    beam.emittance[3] = beam.emittance[3] * beam.beta * beam.total_energy / 2.99792458e8
    return nothing
end

function initilize_6DGaussiandist!(beam::BunchedBeam, optics::AbstractOptics4D, lmap::AbstractLongitudinalMap, cutoff::Float64=4.0)
    # 6D Gaussian distribution

    dist=Truncated(Normal(0.0,1.0),-cutoff,cutoff)

    beam.dist.x .= rand.(dist)
    beam.dist.px .= rand.(dist)
    beam.dist.y .= rand.(dist)
    beam.dist.py .= rand.(dist)
    beam.dist.z .= rand.(dist)
    beam.dist.dp .= rand.(dist)

    get_centroid!(beam)

    beam.dist.x .= beam.dist.x .- beam.centroid[1]
    beam.dist.px .= beam.dist.px .- beam.centroid[2]
    beam.dist.y .= beam.dist.y .- beam.centroid[3]
    beam.dist.py .= beam.dist.py .- beam.centroid[4]
    beam.dist.z .= beam.dist.z .- beam.centroid[5]
    beam.dist.dp .= beam.dist.dp .- beam.centroid[6]

    get_2nd_moment!(beam)
    
    eigval,eigvec=eigen(beam.moment2nd)
    
    mscale = eigvec * diagm(1.0 ./ sqrt.(eigval)) * eigvec'
    # # #beam.dist = mscale*beam.dist
    beam.temp1  .= mscale[1,1] .* beam.dist.x .+ mscale[1,2] .* beam.dist.px .+ mscale[1,3] .* beam.dist.y .+ mscale[1,4] .* beam.dist.py .+ mscale[1,5] .* beam.dist.z .+ mscale[1,6] .* beam.dist.dp
    beam.temp2  .= mscale[2,1] .* beam.dist.x .+ mscale[2,2] .* beam.dist.px .+ mscale[2,3] .* beam.dist.y .+ mscale[2,4] .* beam.dist.py .+ mscale[2,5] .* beam.dist.z .+ mscale[2,6] .* beam.dist.dp
    beam.temp3  .= mscale[3,1] .* beam.dist.x .+ mscale[3,2] .* beam.dist.px .+ mscale[3,3] .* beam.dist.y .+ mscale[3,4] .* beam.dist.py .+ mscale[3,5] .* beam.dist.z .+ mscale[3,6] .* beam.dist.dp
    beam.temp4  .= mscale[4,1] .* beam.dist.x .+ mscale[4,2] .* beam.dist.px .+ mscale[4,3] .* beam.dist.y .+ mscale[4,4] .* beam.dist.py .+ mscale[4,5] .* beam.dist.z .+ mscale[4,6] .* beam.dist.dp
    beam.temp5  .= mscale[5,1] .* beam.dist.x .+ mscale[5,2] .* beam.dist.px .+ mscale[5,3] .* beam.dist.y .+ mscale[5,4] .* beam.dist.py .+ mscale[5,5] .* beam.dist.z .+ mscale[5,6] .* beam.dist.dp
    beam.dist.dp .= mscale[6,1] .* beam.dist.x .+ mscale[6,2] .* beam.dist.px .+ mscale[6,3] .* beam.dist.y .+ mscale[6,4] .* beam.dist.py .+ mscale[6,5] .* beam.dist.z .+ mscale[6,6] .* beam.dist.dp
    
    
    beam.dist.x .= beam.temp1 .* sqrt(beam.emittance[1]*optics.optics_x.beta)
    beam.dist.px .= beam.temp2 .* sqrt(beam.emittance[1]/optics.optics_x.beta)
    beam.dist.y .= beam.temp3 .* sqrt(beam.emittance[2]*optics.optics_y.beta)
    beam.dist.py .= beam.temp4 .* sqrt(beam.emittance[2]/optics.optics_y.beta)
    beam.dist.px .+= beam.dist.x .* (optics.optics_x.alpha/optics.optics_x.beta)
    beam.dist.py .+= beam.dist.y .* (optics.optics_y.alpha/optics.optics_y.beta)

    
    #generate longitudinal distribution based on small amplitude approximation
    eta_p=lmap.αc-1.0/beam.gamma^2
    #println(eta_p)
    
    Qs=sqrt(lmap.RF.v*lmap.RF.h*abs(eta_p*cos(lmap.RF.ϕs))/2/π/beam.beta^2/beam.total_energy)
    println("Qs is ", Qs)
    

    #delta_max = sqrt((lmap.RF.v/eta_p/beam.beta^2/beam.total_energy/lmap.RF.h/abs(eta_p))*(-2*cos(lmap.RF.ϕs)+sin(π - 2*lmap.RF.ϕs)))
    #println(delta_max)
 
    emit_deltap_z = beam.emittance[3]*2.99792458e8/beam.beta/beam.total_energy
    invbeta_deltap_z = Qs*lmap.RF.k/lmap.RF.h/abs(eta_p)
    beam.dist.z .= beam.temp5 .* sqrt(emit_deltap_z/invbeta_deltap_z)
    beam.dist.dp .= beam.dist.dp .* sqrt(emit_deltap_z*invbeta_deltap_z)

    #check separatrix
    #beam phase
    beam.temp1 .= beam.dist.z.*(-lmap.RF.k) .+ lmap.RF.ϕs #beam phase: ϕ
    
    #-cos(ϕs)-cos(ϕ)+sin(ϕs)(π-ϕs-ϕ)
    beam.temp2 .= (-cos(lmap.RF.ϕs)) .- cos.(beam.temp1) + sin(lmap.RF.ϕs).* (π - lmap.RF.ϕs .- beam.temp1)
    
    #delta separatrix * factor=0.95, height bucket or max amplitude of momentum deviation
    beam.temp3 .= sqrt.((lmap.RF.v/lmap.RF.h/(eta_p)/π/beam.beta^2/beam.total_energy) .* beam.temp2) .*0.95
    
    #reset the momentum deviation if it is too large
    beam.dist.dp[(abs.(beam.dist.dp) .> (beam.temp3))] .= beam.temp3[(abs.(beam.dist.dp) .> beam.temp3)] .* (rand.() .- 0.5) .*2.0 


    return nothing
end

function histogram1DinZ!(ps6dcoor::AbstractVector{ps6d{T}}, nbins::Int64, inzindex, zhist, zhist_edges) where T
    # histogram in z
    num_macro=length(ps6dcoor.z)
    zmax=maximum(ps6dcoor.z)
    zmin=minimum(ps6dcoor.z)
    zhist .= 0.0
    zhist_edges .= collect(range(zmin-(zmax-zmin)/nbins, zmax+(zmax-zmin)/nbins, length=nbins+1))
    zsep=(zhist_edges[end]-zhist_edges[1])/nbins
    @inbounds for i in 1:num_macro
        ibin = (ps6dcoor.z[i] - zhist_edges[1])/zsep  # number of bin from 0
        dx = round(ibin) - ibin
        binnum=Int64(floor(ibin)+1)
        inzindex[i] = binnum 
        neighbor = binnum + Int64(sign(dx))
        ratio = (0.5 - abs(dx))/0.5
        weight_neighbor = 0.5 * ratio^2
        zhist[binnum] += 1.0 - weight_neighbor
        zhist[neighbor] += weight_neighbor
    end
    return nothing
end


function histogram1DinZ!(beam::BunchedBeam)
    # histogram in z
    histogram1DinZ!(beam.dist, beam.znbin, beam.inzindex, beam.zhist, beam.zhist_edges)
    return nothing
end


        

