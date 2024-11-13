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

turns = 10000
num_particles = 100000 #num_macro

# Create the proton beam at IP

begin
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
    #distz=open(deserialize, "ebeam_dist_z_wake3x.jls")
    
    #initilize_zslice!(estrong, distz, :evennpar, 11)
    #crab_crossing_setup!(estrong, 12.5e-3)
    
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)
    ilb = InvLorentzBoost(12.5e-3)
    #track!(ebeam, ilb)
    ezcrab2 = easyCrabCavity(394.0e6, 12.5e-3, π) #wrong, crabcrossing setup in new EPIC update

    #define SC_lens for SC_kick BE
    
    nSC = 10 #4
    ds = 3800/(nSC+1)   #0.5     #3833.8451/(nSC+1)
    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0) #random number
    #opSC = optics4DUC(10.0, 0.0, 10.0, 0.0)
    phi_advx = LinRange(0, 29.228, nSC+1)
    phi_advy = LinRange(0, 30.210, nSC+1)
    sc_BE = SC_lens(opSC, ds)
    
    TM1 = Array{TransferMap4D, 1}(undef, nSC)
    TM1_inv = Array{Inverse_TransferMap4D, 1}(undef, nSC)
    

    #records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)
    delta_py = zeros(Float64, num_particles)

    records = StructArray{record}(undef, turns)

    #define SC_kick sympl_gridless
    dt = 3800/5
    sc = SPACECHARGE(dt, 10, 10, 0.001, 0.001)

    #beam beam
    lumi=0.0

end


@benchmark begin

    for i in 1:turns
        track!(pbeam, oneturn)
        track!(pbeam, lmap)
        track!(pbeam, pcrabu)
        track!(pbeam, pcrabu2nd)
        track!(pbeam, lb)
        lumi=track!(pbeam, estrong)
        track!(pbeam, invlb)
        track!(pbeam, pcrabd)
        track!(pbeam, pcrabd2nd)
        

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
            
            SC_kick!(sc_BE, pbeam)
            
            TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
            track!(pbeam, TM1_inv[j])
        end
        
        

        #SC Symplectic gridless model
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)


        particles_turns_SC[:, 1, i] .= pbeam.dist.x
        particles_turns_SC[:, 2, i] .= pbeam.dist.px
        particles_turns_SC[:, 3, i] .= pbeam.dist.y
        particles_turns_SC[:, 4, i] .= pbeam.dist.py
        particles_turns_SC[:, 5, i] .= pbeam.dist.z
        particles_turns_SC[:, 6, i] .= pbeam.dist.dp

        delta_px .=  particles_turns[:,2,i] .- particles_turns_SC[:,2,i]
        delta_py .=  particles_turns[:,4,i] .- particles_turns_SC[:,4,i]
    
        get_emittance!(pbeam)
        records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                            pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                            pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)

    end

end


#SC kick plots
plot(particles_turns[:, 1, turns], delta_px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "Δpx")
plot(particles_turns[:, 3, turns], delta_py, seriestype=:scatter, markersize=0.5, xlabel="y", ylabel= "Δpy")
plot(particles_turns[:, 5, turns], delta_px, seriestype=:scatter, markersize=0.5, xlabel="z", ylabel= "Δpx")

#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1,xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px")
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py")

#check parameters
particles_turns
particles_turns_SC
delta_px

######### tune diagram #########
using PyCall
@pyimport nafflib

half=turns ÷ 2

tunex_1=nafflib.multiparticle_tunes(particles_turns_SC[:, 1, 1:half ])
tunex_2=nafflib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=nafflib.multiparticle_tunes(particles_turns_SC[:, 3, 1:half])
tuney_2=nafflib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .5,  color = :jet, clim=(-10,-2), 
aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300)

    #size=(600,500)#, xlim=(0.228, 0.243), ylim=(0.21, 0.225), xlabel="Horizontal Tune", ylabel="Vertical Tune", label="something", dpi=300)


maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)
##### FFT #####
using LaTeXStrings

plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns")
plot(records.sx, label = L"$σ_x$", xlabel="# Turns")
plot(records.ex, label = L"$ε_x$", xlabel="# Turns")
plot(records.sz, label = L"$σ_z$", xlabel="# Turns")

#remember to modify with the number of turns!
fftx=fft(@view records.cx[5001:end]) # N/2 + 1 to N/2
ffty=fft(@view records.cy[5001:end])
fftz=fft(@view records.cz[5001:end])
fftf=fftfreq(5000, 1.0)

turns = 10000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0,0.25), legendfontsize = 20, label = L"$σ_x$")
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]), xlim=(0,0.25), legendfontsize = 20, label = L"$σ_y$")

maximum(real(fftx))

#save data
open("pbeam_records_tracking_10kturns_100kpart_BB.jls", "w") do f
    serialize(f, records)
end

open("pbeam_records_tracking_100kturns_100kpart_BB_dist.jls", "w") do f
    serialize(f, particles_turns_SC)
end



records = open(deserialize, "pbeam_records_tracking_10kturns_100kpart_BB.jls")
particles_turns_SC = open(deserialize, "pbeam_records_tracking_10kturns_100kpart_SCBE_BB_dist.jls")


# function to save in file .txt
function save_to_file(filename, vector)
    open(filename, "w") do file
        for elem in vector
            write(file, "$elem\n")
        end
    end
end

filename = "test_records.txt"
save_to_file(filename, records)


#calculate average:
#avg_emix = mean(records_1x.ex)
#avg_emiy = mean(records_1y.ey)
avg_sx = mean(records_1x.sx)
avg_sy = mean(records_1y.sy)
avg_sz = mean(records_1y.sz)

factor_tune = 2*pbeam.particle.classrad0*pbeam.num_particle/pbeam.beta^2/pbeam.gamma^3/(2*pi)^1.5/avg_sz
Δνx = factor_tune*3833.8451 /(avg_sx*(avg_sx+avg_sy))*10.0 #last number: beta opSC
Δνy = factor_tune*3833.8451 /(avg_sy*(avg_sx+avg_sy))*10.0

avg_emix = mean(records_1x.ex)
avg_emiy = mean(records_1y.ey)

#solution maybe for tune shift???
Δνx = factor_tune*3833.8451*sqrt(10)/sqrt(avg_emix)/(sqrt(avg_emix)+sqrt(avg_emiy))
Δνy = factor_tune*3833.8451*sqrt(10)/sqrt(avg_emiy)/(sqrt(avg_emix)+sqrt(avg_emiy))

sqrt(avg_emiy)/sqrt(avg_emix)
Δνx/Δνy

# Save the results into file
#=
open("pbeam_wake3x_tracking.jls", "w") do f
    serialize(f, records)
end

open("pbeam_wake3x_status.jls", "w") do f
    serialize(f, pbeam)
end
=#

# Read file and process it later
#=
records_1x=open(deserialize, "pbeam_wake1x_tracking.jls")
plot(records_1x.sy, label="cz")
records_3x=open(deserialize, "pbeam_wake3x_tracking.jls")
plot!(records_3x.sy, label="cz")
=#
