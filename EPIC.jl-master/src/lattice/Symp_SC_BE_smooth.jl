mutable struct SC_lens <: AbstractSpaceCharge
    optics::AbstractOptics4D
    ds::Float64
    nSC::Int64
    turns::Int64
    sigma_x_SC::Matrix{Float64}
    sigma_y_SC::Matrix{Float64}

    function SC_lens(optics, ds, nSC, turns)
        sigma_x_SC = zeros(Float64, turns, nSC)
        sigma_y_SC = zeros(Float64, turns, nSC)
        new(optics, ds, nSC, turns, sigma_x_SC, sigma_y_SC)
    end
end


function parameters_derivatives(x, y, σx, σy)
    
    sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
    termexp = exp(-x*x/2/σx/σx - y*y/2/σy/σy)
    #term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
    #term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
    
    term1=erfcx((x+1im*y)/sqrtδσ2)
    term2=erfcx((x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
    
    ddx_term1 = 2im/sqrt(pi)/sqrtδσ2 - 2*(x+1im*y)/sqrtδσ2^2*term1
    ddx_term2 = termexp * term2* ((2*x + 2im*y)/sqrtδσ2^2) - termexp* (2im*σy/σx/sqrt(pi)/sqrtδσ2)
    ddx_term = ddx_term1 + ddx_term2

    ddy_term1 =  -2/sqrt(pi)/sqrtδσ2 - 2im*(x+1im*y)/sqrtδσ2^2*term1
    ddy_term2 =  termexp * term2* ((2im*x-2*y)/sqrtδσ2^2) + termexp* (2*σx/σy/sqrt(pi)/sqrtδσ2)
    ddy_term = ddy_term1 + ddy_term2 

    dEx_dx = imag(ddx_term) #dEx_dx = dEy_dy so, imag(ddx_term) = real(ddy_term)
    dEx_dy = imag(ddy_term) #dEx_dy = -dEy_dx so, imag(ddy_term) = real(ddx_term)
    dEy_dx = real(ddx_term)
    dEy_dy = real(ddy_term)

    #println(ddx_term1, "\n", ddx_term2,"\n", ddy_term,"\n", ddy_term2)
    #println("dEx_dx:", dEx_dx, "\n", "dEx_dy:", dEx_dy, "\n", "dEy_dx:", dEy_dx, "\n", "dEy_dy:", dEy_dy, "\n")
    
return dEx_dx, dEx_dy, dEy_dx, dEy_dy 
end


function track!(σx, σy, σz, temp1, bbSC::BunchedBeam, factorSC::Float64)

    # temporary solution
    dpx_dz = zeros(length(bbSC.dist.x))
    dpy_dz = zeros(length(bbSC.dist.x))

    dEx_dx = zeros(length(bbSC.dist.x))
    dEx_dy = zeros(length(bbSC.dist.x))
    dEy_dx = zeros(length(bbSC.dist.x))
    dEy_dy = zeros(length(bbSC.dist.x))


    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]
    @inbounds Threads.@threads :static for j in eachindex(bbSC.dist.x)

    
        dEx_dx[j], dEx_dy[j], dEy_dx[j], dEy_dy[j] = parameters_derivatives(bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)

        # gaussian function for particle distribution, lambda_z
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*bbSC.dist.z[j]^2/σz^2)

        dpx_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][1]* bbSC.dist.z[j]/σz^2
        dpy_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][2]* bbSC.dist.z[j]/σz^2

        # delta p_/p_0
        #bbSC.dist.px[j] += factorSC* temp1[j]* (fieldvec_thread[Threads.threadid()][1]+
                        #dEx_dx[j] * bbSC.dist.x[j] + dEx_dy[j] * bbSC.dist.y[j] + dpx_dz[j] * bbSC.dist.z[j])

        #bbSC.dist.py[j] += factorSC* temp1[j]* (fieldvec_thread[Threads.threadid()][2]+
                        #dEy_dx[j] * bbSC.dist.x[j] + dEy_dy[j] * bbSC.dist.y[j] + dpy_dz[j] * bbSC.dist.z[j])

        bbSC.dist.px[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        bbSC.dist.py[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][2]

        bbSC.dist.dp[j] += dpx_dz[j]*bbSC.dist.x[j] + dpy_dz[j]*bbSC.dist.y[j]

        #=
        println("dpx_dz is:", dpx_dz)
        println("dpy_dz is:", dpy_dz)
        println("dEx_dx:", dEx_dx, "\n", "dEx_dy:", dEx_dy, "\n", "dEy_dx:", dEy_dx, "\n", "dEy_dy:", dEy_dy, "\n")
        println("temp1:", temp1)=#

    end
end



#SC kick for the first turn, usefull for collecting beamsize at each j SC-location
function SC_kick!(SC::SC_lens, bbSC::BunchedBeam)
    
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds
    
    get_emittance!(bbSC)
    σx = bbSC.beamsize[1]
    σy = bbSC.beamsize[3] #sqrt(bbSC.emittance[2] * betay) # beamsize y at SC_point
    σz = bbSC.beamsize[5]

    track!(σx, σy, σz, bbSC.temp1, bbSC, factorSC)
    
    return σx, σy
    
end


#function for next k-turns
function SC_kick!(SC::SC_lens, bbSC::BunchedBeam, w::Float64, sx, sy)
    
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds

    get_emittance!(bbSC)

    σx_NEW = sx * w + bbSC.beamsize[1] * (1-w)
    σy_NEW = sy * w + bbSC.beamsize[3] * (1-w)

    σz = bbSC.beamsize[5] #sqrt(bbSC.emittance[2] * betay) # beamsize y at SC_point

    track!(σx_NEW, σy_NEW, σz, bbSC.temp1, bbSC, factorSC)

    return σx_NEW, σy_NEW
    
end

#=
println("Before track! σx: ", σx_NEW, " σy: ", σy_NEW, " σz: ", σz)
track!(σx_NEW, σy_NEW, σz, bbSC.temp1, bbSC, factorSC)
println("After track! σx: ", σx_NEW, " σy: ", σy_NEW, " σz: ", σz)
=#