struct SC_lens <: AbstractSpaceCharge
    optics::AbstractOptics4D
    ds::Float64   #ring arc length
    SC_lens(optics, ds) = new(optics, ds)
end

function track!(σx, σy, σz, temp1, bbSC::BunchedBeam, factorSC::Float64)

    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]
    @inbounds Threads.@threads :static for j in eachindex(bbSC.dist.x)

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)
        
        # gaussian function for particle distribution, lambda_z  (bbSC.dist.z[j]= temp1[j])
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*bbSC.dist.z[j]^2/σz^2)

        # delta p_/p_0
        bbSC.dist.px[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        bbSC.dist.py[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][2]

    end
end


function SC_kick!(SC::SC_lens, bbSC::BunchedBeam)
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds
    
    get_emittance!(bbSC)
    betax = SC.optics.optics_x.beta
    betay = SC.optics.optics_y.beta
    σz = bbSC.beamsize[5]
    σx = sqrt(bbSC.emittance[1] * betax)  # beamsize x at SC_point
    σy = sqrt(bbSC.emittance[2] * betay) #beamsize y at SC_point
    
    track!(σx, σy, σz, bbSC.temp1, bbSC, factorSC)

end











#Analytical tune shift
#=emi_factorx = 1.0/bbSC.emittance[1]/(1.0+sqrt(bbSC.emittance[2]*betay/betax/bbSC.emittance[1]))
emi_factory = 1.0/bbSC.emittance[2]/(1.0+sqrt(bbSC.emittance[1]*betax/betay/bbSC.emittance[2]))

Δνx = []
Δνy = []
Δνx = push!(Δνx, factorSC/SC.ds*emi_factorx/(2*π))
Δνy = push!(Δνy, factorSC/SC.ds*emi_factory/(2*π))

#save in .txt
open("tunex_file.txt", "a") do io
    writedlm(io, Δνx, "\n")
end
open("tuney_file.txt", "a") do io
    writedlm(io, Δνy, "\n")
end=#
#println("Tune shift x,y:", Δνx, " " ,Δνy)

#Equations for tune shift ps
#Δνx = ro*N/beta^2/gamma^3/(2*pi)/Ex * 2/(1+sqrt(betay*Ey/betax/Ex))
#Δνy = ro*N/beta^2/gamma^3/(2*pi)/Ey * 2/(1+sqrt(betax*Ex/betay/Ey))
