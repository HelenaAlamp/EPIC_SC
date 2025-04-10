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
    

end


# Read the data from JuTrack
begin

    data_fodo= readdlm("test/data_fodocell_Ji_cut_in2.txt") #dmux and dmuy have been div by /(2pi)

    betax = round.(data_fodo[1, :], digits=6)
    betay = round.(data_fodo[2,:], digits=6)
    alphax = round.(data_fodo[3,:], digits=6)
    alphay = round.(data_fodo[4,:], digits=6)
    dmux = round.(data_fodo[5,:], digits=6)
    dmuy = round.(data_fodo[6,:], digits=6)

    emix = round((0.0006550846646644346^2)/ betax[10], digits=15)
    emiy = round((0.0007168508938368251^2)/ betay[10], digits=15) #0.0006566409601378707 #0.0007168508938368251 (ex=ey)
    gammax = (1+alphax[10]^2)/betax[10]
    gammay = (1+alphay[10]^2)/betay[10]
    s_px = sqrt(emix*gammax)
    s_py = sqrt(emiy*gammay)

end

#sqrt(emix*betay[5]) to have same emmittance
emix
emiy

emix/emiy

#plot FODO cell parameters
s = [0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.75, 0.8, 0.9, 1]
plot(s, [betax, betay], xticks = s, xlabel="s (m)", linewidth=2, label=["β_x" "β_y"])
plot(s, [alphax, alphay], xticks = s, xlabel="s (m)", linewidth=2, label=["α_x" "α_y"])
plot(s, [dmux, dmuy], xticks = s, xlabel="s (m)", linewidth=2, label=["μ_x" "μ_y"])


#0.25e11 numparticles to have same ji s perveance (e-7)
begin
    turns = 10000
    num_particles = 5000
    
    # generate the bunch
    pbeam = BunchedBeam(PROTON, 0.25e10, 1e9, num_particles, [emix, emiy , 2.3e-5])
    opbeam = optics4DUC(betax[10], alphax[10], betay[10], alphay[10]) #calcolato beta da sigma e emi? - corrisponde a Q3
    mainRF = AccelCavity(100e6, 85e3, 1, 0.0)
    αc=0.0
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opbeam, lmap)
    oneturn = TransferMap4DChrom(opbeam, dmux[10], dmuy[10], 0.0, 0.0)

    # single element optics
    opD11 = optics4DUC(betax[1], alphax[1], betay[1], alphay[1])
    opD12 = optics4DUC(betax[2], alphax[2], betay[2], alphay[2])
    opQ11 = optics4DUC(betax[3], alphax[3], betay[3], alphay[3])
    opQ12 = optics4DUC(betax[4], alphax[4], betay[4], alphay[4])
    opD21 = optics4DUC(betax[5], alphax[5], betay[5], alphay[5])
    opD22 = optics4DUC(betax[6], alphax[6], betay[6], alphay[6])
    opQ21 = optics4DUC(betax[7], alphax[7], betay[7], alphay[7])
    opQ22 = optics4DUC(betax[8], alphax[8], betay[8], alphay[8])
    opD31 = optics4DUC(betax[9], alphax[9], betay[9], alphay[9])
    opD32 = optics4DUC(betax[10], alphax[10], betay[10], alphay[10])

    
    # single element phase advance
    Qx_D11 = dmux[1]
    Qx_D12 = dmux[2]
    Qx_Q11 = dmux[3]
    Qx_Q12 = dmux[4]
    Qx_D21 = dmux[5]
    Qx_D22 = dmux[6]
    Qx_Q21 = dmux[7]
    Qx_Q22 = dmux[8]
    Qx_D31 = dmux[9]
    Qx_D32 = dmux[10]

    Qy_D11 = dmuy[1]
    Qy_D12 = dmuy[2]
    Qy_Q11 = dmuy[3]
    Qy_Q12 = dmuy[4]
    Qy_D21 = dmuy[5]
    Qy_D22 = dmuy[6]
    Qy_Q21 = dmuy[7]
    Qy_Q22 = dmuy[8]
    Qy_D31 = dmuy[9]
    Qy_D32 = dmuy[10]

    #define SC_lens for SC_kick BE
    l_D1 = 0.2
    l_Q1 = 0.1
    l_D2 = 0.4
    l_Q2 = 0.1
    l_D3 = 0.2

    c = 2
    nsc_D1 = c #4
    nsc_Q1 = c #20
    nsc_D2 = c #8
    nsc_Q2 = c #20
    nsc_D3 = c #4

    ds_D1 = l_D1/nsc_D1
    ds_Q1 = l_Q1/nsc_Q1
    ds_D2 = l_D2/nsc_D2
    ds_Q2 = l_Q2/nsc_Q2
    ds_D3 = l_D3/nsc_D3

    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0) # non mi serve in single element

    sc_D1 = SC_lens(opSC, ds_D1, nsc_D1, turns)
    sc_Q1 = SC_lens(opSC, ds_Q1, nsc_Q1, turns)
    sc_D2 = SC_lens(opSC, ds_D2, nsc_D2, turns)
    sc_Q2 = SC_lens(opSC, ds_Q2, nsc_Q2, turns)
    sc_D3 = SC_lens(opSC, ds_D3, nsc_D3, turns)

    #records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    
    delta_px = zeros(Float64, num_particles)
    delta_py = zeros(Float64, num_particles)
    delta_pz = zeros(Float64, num_particles)

    records = StructArray{record}(undef, turns)

    SC = true

    sx = zeros(Float64, turns, 5)
    ex = zeros(Float64, turns, 5)
    sy = zeros(Float64, turns, 5)
    
end


Qx_D31
Qx_D32

get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize

#unbunched beam SC tune shift
ΔQx_SC_emi = -pbeam.particle.classrad0*pbeam.num_particle/(2*pi)/pbeam.gamma^3/pbeam.beta^2/pbeam.emittance[1]
ΔQy_SC_emi = -pbeam.particle.classrad0*pbeam.num_particle/(2*pi)/pbeam.gamma^3/pbeam.beta^2/pbeam.emittance[2]

#bunched beam SC tune shift
ds = 1
sigma_z = 1 #unbunched beam
ΔQx_SC = 2*pbeam.particle.classrad0*pbeam.num_particle/(2*pi)^(3/2)/pbeam.gamma^3/pbeam.beta^2* 
            betax[5]/pbeam.beamsize[1]/(pbeam.beamsize[1]+pbeam.beamsize[3])*ds/sigma_z
ΔQy_SC = 2*pbeam.particle.classrad0*pbeam.num_particle/(2*pi)^(3/2)/pbeam.gamma^3/pbeam.beta^2* 
            betay[5]/pbeam.beamsize[3]/(pbeam.beamsize[1]+pbeam.beamsize[3])*ds/sigma_z


##### compare tune shift ###
#ΔQx-y_SC= 0.29/0.35

#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px")
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py")

begin
    println("Start tracking")
    prog = Progress(turns)

    for i in 1:turns
        
        """particles_turns[:, 1, i] .= pbeam.dist.x
        particles_turns[:, 2, i] .= pbeam.dist.px
        particles_turns[:, 3, i] .= pbeam.dist.y
        particles_turns[:, 4, i] .= pbeam.dist.py
        particles_turns[:, 5, i] .= pbeam.dist.z
        particles_turns[:, 6, i] .= pbeam.dist.dp """
        
        ###### SC tracking - keep same reference - with invTM #####
        if SC == true

            sc_in_element(pbeam, opbeam, opD11, sc_D1, Qx_D11, Qy_D11, "N", "1")
            sc_in_element(pbeam, opbeam, opD12, sc_D1, Qx_D12, Qy_D12, "N", "1")
            """get_emittance!(pbeam)
            sx[i, 1] = pbeam.beamsize[1]
            ex[i, 1] = pbeam.emittance[1] 
            sy[i, 1] = pbeam.beamsize[3]"""

            sc_in_element(pbeam, opbeam, opQ11, sc_Q1, Qx_Q11, Qy_Q11, "N", "1")
            sc_in_element(pbeam, opbeam, opQ12, sc_Q1, Qx_Q12, Qy_Q12, "N", "1")
            """get_emittance!(pbeam)
            sx[i, 2] = pbeam.beamsize[1]
            ex[i, 2] = pbeam.emittance[1] 
            sy[i, 2] = pbeam.beamsize[3]"""

            sc_in_element(pbeam, opbeam, opD21, sc_D2, Qx_D21, Qy_D21, "N", "1")
            sc_in_element(pbeam, opbeam, opD22, sc_D2, Qx_D22, Qy_D22, "N", "1")
            """get_emittance!(pbeam)
            sx[i, 3] = pbeam.beamsize[1]
            ex[i, 3] = pbeam.emittance[1] 
            sy[i, 3] = pbeam.beamsize[3]"""

            sc_in_element(pbeam, opbeam, opQ21, sc_Q2, Qx_Q21, Qy_Q21, "N", "1")
            sc_in_element(pbeam, opbeam, opQ22, sc_Q2, Qx_Q22, Qy_Q22, "N", "1")
            """get_emittance!(pbeam)
            sx[i, 4] = pbeam.beamsize[1]
            ex[i, 4] = pbeam.emittance[1] 
            sy[i, 4] = pbeam.beamsize[3]"""


            sc_in_element(pbeam, opbeam, opD31, sc_D3, Qx_D31, Qy_D31, "N", "1")
            sc_in_element(pbeam, opbeam, opD32, sc_D3, Qx_D32, Qy_D32, "Y", "1")
            """get_emittance!(pbeam)
            sx[i, 5] = pbeam.beamsize[1]
            ex[i, 5] = pbeam.emittance[1] 
            sy[i, 5] = pbeam.beamsize[3]"""
        end
        
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
                            pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3])

        next!(prog)
    end
end

#beam distribution plots
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1, xlabel="x", ylabel= "px", ylim = (-0.005, 0.005), xlim = (-0.002, 0.002))
plot(pbeam.dist.y, pbeam.dist.py, seriestype=:scatter, markersize=1, xlabel="y", ylabel= "py", ylim = (-0.003, 0.004), xlim = (-0.003, 0.002))

# plot some particles last turns
pp = rand(1:1000, 5)

plot(particles_turns_SC[pp, 1, 8000:10000]', particles_turns_SC[pp, 2, 8000:10000]', seriestype=:scatter, markersize=2, xlabel="x", ylabel= "px", labels=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5"])
plot(particles_turns_SC[pp, 3, 8000:10000]', particles_turns_SC[pp, 4, 8000:10000]', seriestype=:scatter, markersize=2, xlabel="y", ylabel= "py", labels=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5"])


# SC kicks phase space
maximum(pbeam.dist.dp)
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1, xlabel="z", ylabel= "δ")
plot(pbeam.dist.x, delta_px, seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "Δpx")
plot(pbeam.dist.y, delta_py, seriestype=:scatter, markersize=0.9, xlabel="y", ylabel= "Δpy")

#beam size - beta function check
plot(sx, ylabel = "σ_x", xlabel= "turns", legendfontsize = 10)
plot(sy, ylabel = "σ_y", xlabel= "turns", legendfontsize = 10)
calc_betax = sx.^2 ./ex
plot(calc_betax[1:1000, :], ylabel = "β_x", xlabel= "turns", legendfontsize = 10)

emix_growth = (records.ex[10000]-records.ex[1])/records.ex[1]*100
betax_var = (mean(calc_betax[:,4]) - betax[5])/betax[5]*100
emiy_growth = (records.ey[10000]-records.ey[1])/records.ey[1]*100


plot(records.ex, label = L"$ε_x$", xlabel="# Turns", legendfontsize = 15)
plot(records.ey, label = L"$ε_y$", xlabel="# Turns", legendfontsize = 15)
plot(records.sx, label = L"$σ_x$", xlabel="# Turns", legendfontsize = 15)
plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns", legendfontsize = 15)
plot(records.sz, label = L"$σ_z$", xlabel="# Turns", legendfontsize = 15)

# multiparticle tune
using PyCall
@pyimport NAFFlib

half=turns ÷ 2
half = 9000

tunex_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, 8000:half])
tunex_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, 8000:half])
tuney_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .9,  color = :jet, clim=(-10,-2), 
    aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 2", dpi=300)
    #size=(600,500)#, xlim=(0.228, 0.243), ylim=(0.21, 0.225), xlabel="Horizontal Tune", ylabel="Vertical Tune", label="something", dpi=300)


maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)


##### FFT #####
using LaTeXStrings

fftx=fft(@view records.cx[8001:9000]) # N/2 + 1 to N/2
fftx1=fft(@view records.cx[9001:10000])

ffty=fft(@view records.cy[8001:9000])
ffty1=fft(@view records.cy[9001:10000])

# separated tunes last 2k
fftf=fftfreq(1000, 1.0)
plot(@view(fftf[1:1000]), abs.(@view fftx[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$x$ - I group", xlim=(0.1,0.35))
plot!(@view(fftf[1:1000]), abs.(@view fftx1[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$x$ - II group")

plot(@view(fftf[1:1000]), abs.(@view ffty[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$y$ - I group", xlim=(0.1,0.35))
plot!(@view(fftf[1:1000]), abs.(@view ffty1[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$y$ - II group")

maximum(abs.(@view fftx[1:1000]))
maximum(abs.(@view fftx1[1:1000]))
maximum(abs.(@view ffty[1:1000]))
maximum(abs.(@view ffty1[1:1000]))


#data = hcat(fftx, fftx1, ffty, ffty1, fftf)
#writedlm("fodocell_FFT_last2k_check_peaks.txt", data)

#fftz=fft(@view records.cz[501:end])
# total half turns tune
fftf=fftfreq(5000, 1.0)
turns = 10000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]),  legendfontsize = 15, xlabel=L"tune", label = L"$centroid - x$", xlim=(0.1,0.35)) #xlim=(0.02,0.035),
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]),  legendfontsize = 15, xlabel=L"tune", label = L"$centroid - y$", xlim=(0.1,0.35)) #xlim=(0.01,0.025),
plot(@view(fftf[1:turns÷2]), abs.(@view fftz[1:turns÷2]),  legendfontsize = 15, xlabel=L"tune", label = L"$centroid - y$") #xlim=(0.01,0.025),


sx_filtered = records.sx.-mean(records.sx)
sy_filtered = records.sy.-mean(records.sy)

fftx=fft(@view sx_filtered[501:end]) # N/2 + 1 to N/2
ffty=fft(@view sy_filtered[501:end])
fftf=fftfreq(500, 1.0)



### peaks from NAFFlib ###
#Q, A = NAFFlib.get_tunes(x, N, order, interpolation)
Qx, Ax = NAFFlib.get_tunes(particles_turns_SC[1000, 1, 8000:9000], 4)
Qx1, Ax1 = NAFFlib.get_tunes(particles_turns_SC[1000, 1, 9001:10000], 4)

scatter(Qx, abs.(Ax), xlabel="Horizontal Tune", ylabel= "Amplitude", label = L"$x$ - I group", legendfontsize = 10)
scatter!(Qx1, abs.(Ax1), label = L"$x$ - II group", legendfontsize = 10)

Qy, Ay = NAFFlib.get_tunes(particles_turns_SC[1000, 3, 8000:9000], 4)
Qy1, Ay1 = NAFFlib.get_tunes(particles_turns_SC[1000, 3, 9001:10000], 4)

scatter(Qy, abs.(Ay), xlabel="Vertical Tune", ylabel= "Amplitude", label = L"$y$ - I group", legendfontsize = 10)
scatter!(Qy1, abs.(Ay1), label = L"$y$ - II group", legendfontsize = 10)



######## fft single particle #######
fftx=fft(@view particles_turns_SC[500, 1, 8001:9000]) # N/2 + 1 to N/2
fftx1=fft(@view particles_turns_SC[500, 1, 9001:10000])

ffty=fft(@view particles_turns_SC[500, 3, 8001:9000])
ffty1=fft(@view particles_turns_SC[500, 3, 9001:10000])

# separated tunes last 2k
fftf=fftfreq(1000, 1.0)
plot(@view(fftf[1:1000]), abs.(@view fftx[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$x$ - I group", xlim=(0.15,0.35))
plot!(@view(fftf[1:1000]), abs.(@view fftx1[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$x$ - II group")

plot(@view(fftf[1:1000]), abs.(@view ffty[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$y$ - I group", xlim=(0.1,0.35))
plot!(@view(fftf[1:1000]), abs.(@view ffty1[1:1000]), legendfontsize = 15, xlabel=L"tune", label = L"$y$ - II group")

maximum(abs.(@view fftx[1:1000]))
maximum(abs.(@view fftx1[1:1000]))
maximum(abs.(@view ffty[1:1000]))
maximum(abs.(@view ffty1[1:1000]))



peak_amplitude = maximum(abs.(ffty1))

peak_index = findfirst(x -> abs(x) == peak_amplitude, abs.(ffty1))
peak_frequency = fftf[peak_index]

println("Peak amplitude: ", peak_amplitude)
println("Peak index: ", peak_index)
println("Peak frequency (tune): ", peak_frequency)