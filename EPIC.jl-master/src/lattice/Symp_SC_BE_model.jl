
struct SC_lens <: AbstractSpaceCharge
    optics::AbstractOptics4D
    ds::Float64   #ring arc length
    SC_lens(optics, ds) = new(optics, ds)
end


function parameters_derivatives(x, y, σx, σy)

    sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
    term3 = Complex(exp(-(x+im*y)^2 /sqrtδσ2^2))
    term4 = Complex(exp(-(x*σy/σx + 1im*y*σx/σy)^2/sqrtδσ2^2))
    term2 = erfcx((x*σy/σx+1im*y*σx/σy)/sqrtδσ2)    #div term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
    termexp = exp(-x*x/2/σx/σx-y*y/2/σy/σy)

    ddx_term = (2/sqrtδσ2/sqrt(pi)* term3) + termexp* ((x/σx^2)* term2 - 
        term4* (2*σy/σx/sqrt(pi)/sqrtδσ2))

    ddy_term = (2im/sqrtδσ2/sqrt(pi) * term3) - termexp*((y/σy^2)* term2 +
        term4* (2im*σx/σy/sqrt(pi)/sqrtδσ2))


    dEx_dx = imag(ddx_term)
    dEy_dx = real(ddx_term)
    dEx_dy = imag(ddy_term)
    dEy_dy = real(ddy_term)

return dEx_dx, dEx_dy, dEy_dx, dEy_dy
end


function track!(σx, σy, σz, temp1, bbSC::BunchedBeam, factorSC::Float64)

    dpx_dz = zeros(length(bbSC.dist.x))
    dpy_dz = zeros(length(bbSC.dist.x))

    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]
    @inbounds Threads.@threads :static for j in eachindex(bbSC.dist.x)
    

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)

        # gaussian function for particle distribution, lambda_z
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*bbSC.dist.z[j]^2/σz^2)

        dEx_dx, dEx_dy, dEy_dx, dEy_dy = parameters_derivatives(bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)

        dpx_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][1]* bbSC.dist.z[j]/σz^2
        dpy_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][2]* bbSC.dist.z[j]/σz^2

        # delta p_/p_0
        bbSC.dist.px[j] += factorSC* temp1[j]* (fieldvec_thread[Threads.threadid()][1]+
                        dEx_dx * bbSC.dist.x[j] + dEx_dy * bbSC.dist.y[j] + dpx_dz[j] * bbSC.dist.z[j])

        bbSC.dist.py[j] += factorSC* temp1[j]* (fieldvec_thread[Threads.threadid()][2]+
                        dEy_dx * bbSC.dist.x[j] + dEy_dy * bbSC.dist.y[j] + dpy_dz[j] * bbSC.dist.z[j])

        
        #dp_dx = dpx_dz
        #dp_dx = dpy_dz
        bbSC.dist.dp[j] += dpx_dz[j]*bbSC.dist.x[j] + dpy_dz[j]*bbSC.dist.y[j]

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
