using Revise, EPIC
using Test
using BenchmarkTools
using Plots
using DelimitedFiles
using Serialization
using StructArrays
using FFTW
using LaTeXStrings

theme(:ggplot2)
#theme(:vibrant)

struct record
    cx::Float64
    cy::Float64
    cz::Float64
    sx::Float64
    sy::Float64
    sz::Float64
    ex::Float64
    ey::Float64
    ez::Float64
    lumi::Float64
end


begin
    num_particles = 1000
    turns = 1000

    pbeam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-2])  #3.7e-2 -> 6cm
    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0)
    mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
    αc=1.5e-3
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)

    # define crab cavity
    overcrab=1.0
    crab_ratio=0.333
    pcrabu = easyCrabCavity(197.0e6, overcrab*12.5e-3*(1+crab_ratio))
    pcrabu2nd = easyCrabCavity(197.0e6*2.0, -overcrab*12.5e-3*crab_ratio)
    
    pcrabd = easyCrabCavity(197.0e6, -overcrab*12.5e-3*(1+crab_ratio))
    pcrabd2nd = easyCrabCavity(197.0e6*2.0, overcrab*12.5e-3*crab_ratio)

    # define Lorentz boost
    lb = LorentzBoost(12.5e-3)
    invlb = InvLorentzBoost(12.5e-3)

    # define oneturnmap
    tunex=0.228
    tuney=0.21
    oneturn = TransferMap4DChrom(opIPp, tunex, tuney, 1.0, 1.0)
    opIPe = optics4DUC(0.45, 0.0, 0.056, 0.0)
    estrong = StrongGaussianBeam(ELECTRON, 1.72e11, 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 5)
    
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)
    ilb = InvLorentzBoost(12.5e-3)
    ezcrab2 = easyCrabCavity(394.0e6, 12.5e-3, π) #wrong, crabcrossing setup in new EPIC update

    #define SC_lens for SC_kick BE
    nSC = 5
    ds = 3800 / (nSC + 1)
    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0)
    turns_rec = 200
    sc_BE = SC_lens(opSC, ds, nSC, turns_rec)
    w = 0.9

    phi_advx = LinRange(0, 29.228, nSC+1)
    phi_advy = LinRange(0, 30.210, nSC+1)
    
    TM1 = Array{TransferMap4D, 1}(undef, nSC)
    TM1_inv = Array{Inverse_TransferMap4D, 1}(undef, nSC)
    
    #records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)
    records = StructArray{record}(undef, turns)
    
    lumi=0.0

end

# record beam size under SC effects at each nSC location
for i in 1:turns_rec

    track!(pbeam, oneturn)

    for j in 1:nSC
        
        TM1[j] = TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1[j])
        
        sc_BE.sigma_x_SC[i,j], sc_BE.sigma_y_SC[i,j] = SC_kick!(sc_BE, pbeam)
        
        TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1_inv[j])
    end

end


# record in file previous SC beam size
open("pbeam_records_sx_SC_Sympl_TEST.jls", "w") do f
    serialize(f, sc_BE.sigma_x_SC)
end

open("pbeam_records_sy_SC_Sympl_TEST.jls", "w") do f
    serialize(f, sc_BE.sigma_y_SC)
end

plot(sc_BE.sigma_x_SC, label = L"$σ_x$", xlabel="# Turns", )
plot(sc_BE.sigma_y_SC, label = L"$σ_y$", xlabel="# Turns", )

# average beam size for next turn
for j in 1:nSC
    sc_BE.sigma_x_SC[j] = mean(sc_BE.sigma_x_SC[:,j])
    sc_BE.sigma_y_SC[j] = mean(sc_BE.sigma_y_SC[:,j])
end


############## start tracking  #################

for i in 1:turns
    track!(pbeam, oneturn)
    #=track!(pbeam, lmap)
    track!(pbeam, pcrabu)
    track!(pbeam, pcrabu2nd)
    track!(pbeam, lb)
    lumi=track!(pbeam, estrong)
    track!(pbeam, invlb)
    track!(pbeam, pcrabd)
    track!(pbeam, pcrabd2nd)=#
    
    particles_turns[:, 1, i] .= pbeam.dist.x
    particles_turns[:, 2, i] .= pbeam.dist.px
    particles_turns[:, 3, i] .= pbeam.dist.y
    particles_turns[:, 4, i] .= pbeam.dist.py
    particles_turns[:, 5, i] .= pbeam.dist.z
    particles_turns[:, 6, i] .= pbeam.dist.dp


    #SC BE model 
    for j in 1:(nSC)
        
        TM1[j] = TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1[j])
        
        sc_BE.sigma_x_SC[j], sc_BE.sigma_y_SC[j] = SC_kick!(sc_BE, pbeam, w, sc_BE.sigma_x_SC[j], sc_BE.sigma_y_SC[j])
        

        TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1_inv[j])
    end

    particles_turns_SC[:, 1, i] .= pbeam.dist.x
    particles_turns_SC[:, 2, i] .= pbeam.dist.px
    particles_turns_SC[:, 3, i] .= pbeam.dist.y
    particles_turns_SC[:, 4, i] .= pbeam.dist.py
    particles_turns_SC[:, 5, i] .= pbeam.dist.z
    particles_turns_SC[:, 6, i] .= pbeam.dist.dp

    get_emittance!(pbeam)
    records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                        pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                        pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)

    delta_px .=  particles_turns[:,2,i] .- particles_turns_SC[:,2,i]
    
    
end

plot(pbeam.dist.x, delta_px, seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "Δpx")



#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1,xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px")
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py")


#########################################

plot(records.sx, label = L"$σ_x$", xlabel="# Turns")
plot(records.ex, label = L"$ε_x$", xlabel="# Turns")
plot(records.sz, label = L"$σ_z$", xlabel="# Turns")
plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns")

##########################################
using PyCall
@pyimport NAFFlib

half=turns ÷ 2

#half = 9000
tunex_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, 1:half ])
tunex_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, 1:half])
tuney_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .9,  color = :jet, clim=(-10,-2), 
aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300)#,  xlim=(0.22, 0.23), ylim=(0.185, 0.215))

maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)

####################################################################

fftx=fft(@view records_after.cx[5001:end]) # N/2 + 1 to N/2
ffty=fft(@view records_after.cy[5001:end])
fftz=fft(@view records_after.cz[5001:end])
fftf=fftfreq(5000, 1.0)

turns = 10000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0.2,0.25), legendfontsize = 10, label = L"$x$ $centroid$", xlabel =L"$tune$")
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]), xlim=(0.2,0.25), legendfontsize = 10, label = L"$y$ $centroid$", xlabel =L"$tune$")


records_after = open(deserialize,"pbeam_records_BBSC_10kturns_10kpart_w05.jls")
plot(records_after.sx, legendfontsize = 10, label = L"$σ_x$", xlabel="# Turns")
plot(records_after.ex, legendfontsize = 10, label = L"$ε_x$", xlabel="# Turns")

