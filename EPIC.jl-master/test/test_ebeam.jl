using Revise, EPIC
using Test
using BenchmarkTools
using Plots
theme(:ggplot2)

begin
    ebeam=BunchedBeam(ELECTRON, 1.72e11, 10e9,  1000, [20e-9, 1.3e-9, 1.36e-4])
    opIPe=optics4DUC(0.45,0.0,0.056,0.0)
    vbase=3.42*8.5e6
    ϕs=10
    vact=vbase/cos(ϕs*π/180.0)
    mainRFe=AccelCavity(591e6, vact, 7560.0, π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    tunez=get_synchrotron_tune(ebeam, lmap)
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)

    get_emittance!(ebeam)
    eqbs=1.0.*ebeam.beamsize

    ezcrab1=easyCrabCavity(394.0e6, 12.5e-3)
    #track!(ebeam, ezcrab1)
    lb=LorentzBoost(12.5e-3)
    #track!(ebeam, lb)
    #scatter(ebeam.dist.z, ebeam.dist.x; markersize=0.1, color=:oslo)
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    
    pstrong=StrongGaussianBeam(PROTON, 0.688e11, 275e9,  opIPp, [95e-6, 8.5e-6, 0.06], 11)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 5.0)
    #lumi=track!(ebeam, pstrong) 
    ilb=InvLorentzBoost(12.5e-3)
    #track!(ebeam, ilb)
    ezcrab2=easyCrabCavity(394.0e6, 12.5e-3, π)
    #track!(ebeam, ezcrab2)
    #scatter(ebeam.dist.x, ebeam.dist.px; markersize=0.1, markercolor=:red)
    RLCwake = LongitudinalRLCWake(20e9, 10.0e3, 1.0)


    oneturn=TransferMap4DChrom(opIPe, tunex, tuney, 2.0, 2.0)
    oneturndamp=OneTurnRadiation(4000.0, 4000.0, 2000.0)
    eqbs


    line1=Lattice([ezcrab1, lb])
    line2=Lattice([ilb, ezcrab2, oneturn, lmap, RLCwake])
    longiline=Lattice([lmap, RLCwake])

    #add here SC calculations?

end

eqbs
dpcopy=1.0.*ebeam.dist.dp

track!(ebeam, RLCwake)

scatter(ebeam.dist.z, ebeam.dist.dp.-dpcopy)
turns=500
coor_rec = Array{Float64, 2}(undef, 6, turns)
size_rec = Array{Float64, 2}(undef, 6, turns)
lumis = Vector{Float64}(undef, turns)
for i in 1:turns
    track!(ebeam, line1)
    lumis[i]=track!(ebeam, pstrong)
    track!(ebeam, oneturndamp, eqbs)
    track!(ebeam, line2)
    #track!(ebeam, longiline)
    get_emittance!(ebeam)
    coor_rec[:,i] .= ebeam.centroid
    size_rec[:,i] .= ebeam.beamsize
end

size_rec


using FFTW

plot(@view size_rec[5,:])
plot!(@view coor_rec[5,:])
fftx=fft(@view coor_rec[1,251:end])
ffty=fft(@view coor_rec[3,251:end])
fftz=fft(@view coor_rec[5,251:end])
fftf=fftfreq(500, 1.0)
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0,0.3))

get_emittance!(ebeam)
ebeam.emittance
ebeam.centroid
ebeam.beamsize

eqbs

plot(ebeam.dist.z, ebeam.dist.dp, seriestype=:scatter, markersize=0.5)
plot(ebeam.dist.x, ebeam.dist.z, seriestype=:scatter, markersize=0.5)
#plot(ebeam.dist.y, ebeamdist.py, seriestype=:scatter, markersize=0.5)


