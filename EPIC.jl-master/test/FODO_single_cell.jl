using Revise, EPIC
using Test
using BenchmarkTools
using Plots
using DelimitedFiles
using Serialization
using StructArrays
using FFTW
using LaTeXStrings
using ProgressMeter
using Distributions

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


"""
# symmetric fodo cell
# on JuTrack:
initial beam size
0.0006550846646644346
0.001503740507204768
0.0006566409601378707
0.0015094494868860952

eta_p: -0.8803545114201268
Qs: 0.009977004359703896
delta_max: 0.04281724330933424 with 2.3e-5 seV and AccelCavity(100e6, 85e3, 1, 0.0), 1GeV
"""

# Read the data from JuTrack
begin

    data_fodo= readdlm("test/data_fodocell_Ji_not_symmetric_09.txt")

    betax = data_fodo[1, :]
    betay = data_fodo[2,:]
    alphax = data_fodo[3,:]
    alphay = data_fodo[4,:]
    dmux = data_fodo[5,:]
    dmuy = data_fodo[6,:]

    emix = round((0.0006550846646644346^2)/ betax[5], digits=12)
    emiy = round((0.0006566409601378707^2)/betay[5], digits=12)
end


#0.25e11 to have same ji s perveance (e-7)
begin
    turns = 1000
    num_particles = 1000
    
    # generate the bunch
    pbeam = BunchedBeam(PROTON, 0.25e11, 1e9, num_particles, [emix, emiy, 2.3e-5])
    opbeam = optics4DUC(betax[5], alphax[5], betay[5], alphay[5])
    mainRF = AccelCavity(100e6, 85e3, 1, 0.0) 
    αc=0.0
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opbeam, lmap)
    #oneturn = TransferMap4DChrom(opbeam, 1.51916/2/pi, 1.51916/2/pi, 0.0, 0.0)

    # single element optics
    opD1 = optics4DUC(betax[1], alphax[1], betay[1], alphay[1])
    opQ1 = optics4DUC(betax[2], alphax[2], betay[2], alphay[2])
    opD2 = optics4DUC(betax[3], alphax[3], betay[3], alphay[3])
    opQ2 = optics4DUC(betax[4], alphax[4], betay[4], alphay[4])
    opD3 = optics4DUC(betax[5], alphax[5], betay[5], alphay[5]) # =opbeam
    
    # single element phase advance
    Qx_D1 = dmux[1]
    Qx_Q1 = dmux[2]
    Qx_D2 = dmux[3]
    Qx_Q2 = dmux[4]
    Qx_D3 = dmux[5]

    Qy_D1 = dmuy[1]
    Qy_Q1 = dmuy[2]
    Qy_D2 = dmuy[3]
    Qy_Q2 = dmuy[4]
    Qy_D3 = dmuy[5]

    # single element transfermap from pbeam to element - 
    # defined in function sc_in_element
    """TM_D1= TransferMap4D(opbeam, opD1, Qx_D1, Qy_D1)
    TM_Q1= TransferMap4D(opbeam, opQ1, Qx_Q1, Qy_Q1)
    TM_D2= TransferMap4D(opbeam, opD2, Qx_D2, Qy_D2)
    TM_Q2= TransferMap4D(opbeam, opQ2, Qx_Q2, Qy_Q2)
    TM_D3= TransferMap4D(opbeam, opD3, Qx_D3, Qy_D3)
    """
    TM_D31 = TransferMap4D(opbeam, opD3, Qx_D3, Qy_D3)

    # single element transfermap from previous element to next element
    """TM_D1= TransferMap4D(opbeam, opD1, Qx_D1, Qy_D1)
    TM_Q1= TransferMap4D(opD1, opQ1, Qx_Q1-Qx_D1, Qy_Q1-Qx_D1)
    TM_D2= TransferMap4D(opQ1, opD2, Qx_D2-Qx_Q1, Qy_D2-Qx_Q1)
    TM_Q2= TransferMap4D(opD2, opQ2, Qx_Q2-Qx_D2, Qy_Q2-Qx_D2)
    TM_D3= TransferMap4D(opQ2, opD3, Qx_D3-Qx_Q2, Qy_D3-Qx_Q2) #D3==opbeam
    """

    #define SC_lens for SC_kick BE
    l_D1 = 0.2
    l_Q1 = 0.1
    l_D2 = 0.4
    l_Q2 = 0.1
    l_D3 = 0.2

    nsc_D1 = 1#4
    nsc_Q1 = 1#20
    nsc_D2 = 1#8
    nsc_Q2 = 1#20
    nsc_D3 = 1#4

    """ds_D1 = l_D1/nsc_D1
    ds_Q1 = l_Q1/nsc_Q1
    ds_D2 = l_D2/nsc_D2
    ds_Q2 = l_Q2/nsc_Q2
    ds_D3 = l_D3/nsc_D3"""

    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0) # no need

    sc_D1 = SC_lens(opSC, l_D1, nsc_D1, turns)
    sc_Q1 = SC_lens(opSC, l_Q1, nsc_Q1, turns)
    sc_D2 = SC_lens(opSC, l_D2, nsc_D2, turns)
    sc_Q2 = SC_lens(opSC, l_Q2, nsc_Q2, turns)
    sc_D3 = SC_lens(opSC, l_D3, nsc_D3, turns)
    
    """sc_phix_D1 = LinRange(0, Qx_D1, nsc_D1)
    sc_phix_Q1 = LinRange(Qx_D1, Qx_Q1, nsc_Q1)
    sc_phix_D2 = LinRange(Qx_Q1, Qx_D2, nsc_D2)
    sc_phix_Q2 = LinRange(Qx_D2, Qx_Q2, nsc_Q2)
    sc_phix_D3 = LinRange(Qx_Q2, Qx_D3, nsc_D3)

    sc_phiy_D1 = LinRange(0, Qy_D1, nsc_D1)
    sc_phiy_Q1 = LinRange(Qy_D1, Qy_Q1, nsc_Q1)
    sc_phiy_D2 = LinRange(Qy_Q1, Qy_D2, nsc_D2)
    sc_phiy_Q2 = LinRange(Qy_D2, Qy_Q2, nsc_Q2)
    sc_phiy_D3 = LinRange(Qy_Q2, Qy_D3, nsc_D3)"""

    #records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    
    delta_px = zeros(Float64, num_particles)
    delta_py = zeros(Float64, num_particles)
    delta_pz = zeros(Float64, num_particles)

    records = StructArray{record}(undef, turns)
    lumi = 0.0
    SC = true

end

get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize

#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px")
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py")

begin
    println("Start tracking")
    prog = Progress(turns)

    for i in 1:turns

        particles_turns[:, 1, i] .= pbeam.dist.x
        particles_turns[:, 2, i] .= pbeam.dist.px
        particles_turns[:, 3, i] .= pbeam.dist.y
        particles_turns[:, 4, i] .= pbeam.dist.py
        particles_turns[:, 5, i] .= pbeam.dist.z
        particles_turns[:, 6, i] .= pbeam.dist.dp    
        
        ###### SC tracking - keep same reference - with invTM #####
        
        if SC == true
            sc_in_element(pbeam, opbeam, opD1, sc_D1, Qx_D1, Qy_D1)
            sc_in_element(pbeam, opbeam, opQ1, sc_Q1, Qx_Q1, Qy_Q1)
            sc_in_element(pbeam, opbeam, opD2, sc_D2, Qx_D2, Qy_D2)
            sc_in_element(pbeam, opbeam, opQ2, sc_Q2, Qx_Q2, Qy_Q2)
            
            track!(pbeam, TM_D31) #opD3=pbeam oneturnmap
            SC_kick!(sc_D3, pbeam) 
            
        end
        

        ###### SC tracking - element by element #####
        """
        if SC == true
            track!(pbeam, TM_D1)
            SC_kick!(sc_D1, pbeam)
            track!(pbeam, TM_Q1)
            SC_kick!(sc_Q1, pbeam)
            track!(pbeam, TM_D2)
            SC_kick!(sc_D2, pbeam)
            track!(pbeam, TM_Q2)
            SC_kick!(sc_Q2, pbeam)
            track!(pbeam, TM_D3)
            SC_kick!(sc_D3, pbeam)
            track!(pbeam, oneturn)
        end
        """

    
        particles_turns_SC[:, 1, i] .= pbeam.dist.x
        particles_turns_SC[:, 2, i] .= pbeam.dist.px
        particles_turns_SC[:, 3, i] .= pbeam.dist.y
        particles_turns_SC[:, 4, i] .= pbeam.dist.py
        particles_turns_SC[:, 5, i] .= pbeam.dist.z
        particles_turns_SC[:, 6, i] .= pbeam.dist.dp


        delta_px .=  particles_turns[:,2,i] .- particles_turns_SC[:,2,i]
        delta_py .=  particles_turns[:,4,i] .- particles_turns_SC[:,4,i]
        delta_pz .=  particles_turns[:,6,i] .- particles_turns_SC[:,6,i]

        get_emittance!(pbeam)
        records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                            pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                            pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)

        next!(prog)
    end
end

##### compare tune shift ######
Qx_D3-Qy_D3
ds = 1
sigma_z = 1 #unbunched beam
ΔQy_SC = 2*pbeam.particle.classrad0*pbeam.num_particle/(2*pi)^(3/2)/pbeam.gamma^3/pbeam.beta^2* 
            betay[5]/pbeam.beamsize[2]/(pbeam.beamsize[1]+pbeam.beamsize[2])*ds/sigma_z

ΔQx_SC = 2*pbeam.particle.classrad0*pbeam.num_particle/(2*pi)^(3/2)/pbeam.gamma^3/pbeam.beta^2* 
            betax[5]/pbeam.beamsize[1]/(pbeam.beamsize[1]+pbeam.beamsize[2])*ds/sigma_z



maximum(pbeam.dist.dp)
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, delta_px, seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "Δpx")
plot(pbeam.dist.y, delta_py, seriestype=:scatter, markersize=0.9, xlabel="y", ylabel= "Δpy")

emi_growth = (records.ex[1000]-records.ex[1])/records.ex[1]*100

plot(records.ex, label = L"$ε_x$", xlabel="# Turns", legendfontsize = 15)
plot(records.ey, label = L"$ε_y$", xlabel="# Turns", legendfontsize = 15)

plot(records.sx, label = L"$σ_x$", xlabel="# Turns", legendfontsize = 15)
plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns", legendfontsize = 15)
plot(records.sz, label = L"$σ_z$", xlabel="# Turns", legendfontsize = 15)


get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize


emit_growth[i] = (new_emit1[i+1, 1] / new_emit1[1, 1] * new_emit1[i+1, 2] / new_emit1[1, 2] - 1) * 100


#savefig("emittance_Ji_fodocell_10kpart_10kturns_3D.pdf")
######### tune diagram #########

using PyCall
@pyimport NAFFlib

half=turns ÷ 2
half=9000

tunex_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, 1:half])
tunex_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, 1:half])
tuney_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .9,  color = :jet, clim=(-10,-2), 
aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300)

    #size=(600,500)#, xlim=(0.228, 0.243), ylim=(0.21, 0.225), xlabel="Horizontal Tune", ylabel="Vertical Tune", label="something", dpi=300)


maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)

#savefig("tune_Ji_fodocell_10kpart_10kturns_3D.pdf")

#SC kick plots
plot(particles_turns[:, 1, turns], delta_px, seriestype=:scatter, markersize=.8, xlabel="x", ylabel= "Δpx")
plot(particles_turns[:, 3, turns], delta_py, seriestype=:scatter, markersize=.8, xlabel="y", ylabel= "Δpy")
plot(particles_turns[:, 5, turns], delta_px, seriestype=:scatter, markersize=.8, xlabel="z", ylabel= "Δpx")

#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px")
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py")

#check parameters
particles_turns
particles_turns_SC
delta_px


##### FFT #####
using LaTeXStrings

fftx=fft(@view records.cx[501:end]) # N/2 + 1 to N/2
ffty=fft(@view records.cy[501:end])
fftz=fft(@view records.cz[501:end])
fftf=fftfreq(500, 1.0)

turns = 1000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]),  legendfontsize = 15, xlabel=L"tune", label = L"$centroid_x$") #xlim=(0.02,0.035),
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]),  legendfontsize = 15, xlabel=L"tune", label = L"$centroid_y$") #xlim=(0.01,0.025),


sx_filtered = records.sx.-mean(records.sx)
sy_filtered = records.sy.-mean(records.sy)

fftx=fft(@view sx_filtered[5001:end]) # N/2 + 1 to N/2
ffty=fft(@view sy_filtered[5001:end])
fftf=fftfreq(5000, 1.0)


turns = 10000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0.4,0.5),  legendfontsize = 15, xlabel=L"tune", label = L"$σ_x$") #xlim=(0.02,0.035),
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]), xlim=(0.4,0.5),  legendfontsize = 15, xlabel=L"tune", label = L"$σ_y$") #xlim=(0.01,0.025),

