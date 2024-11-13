using Revise, EPIC
using Test
using BenchmarkTools
using Plots
theme(:ggplot2)

begin
    pbeam=BunchedBeam(PROTON, 0.688e11, 275e9,  10000, [11.3e-9, 1.0e-9, 0.0374])
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    vbase=15.87e6 
    ϕs=0
    vact=vbase/cos(ϕs*π/180.0)
    mainRFp=AccelCavity(591e6, vact, 7560.0, π-ϕs*π/180.0) #pi-phi_s = - phi_s
    tunex, tuney = 29.228, 30.210
    αc=1.44e-3
    lmap=LongitudinalRFMap(αc, mainRFp)
    tunez=get_synchrotron_tune(pbeam, lmap)
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)

    get_emittance!(pbeam)
    pqbs=1.0.*pbeam.beamsize


    ezcrab1=easyCrabCavity(394.0e6/2, 12.5e-3)
    #track!(pbeam, ezcrab1)
    lb=LorentzBoost(12.5e-3)
    #track!(pbeam, lb)
    #scatter(pbeam.dist.z, pbeam.dist.x; markersize=0.1, color=:oslo)
   
    opIPe=optics4DUC(0.45, 0.0, 0.056, 0.0)  #add here from above
    estrong=StrongGaussianBeam(ELECTRON, 1.72e11, 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 5) #< nzslice
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)
    #lumi=track!(pbeam, estrong) 
    ilb=InvLorentzBoost(12.5e-3)
    #track!(pbeam, ilb)
    ezcrab2=easyCrabCavity(394.0e6/2, 12.5e-3, π)
    #track!(pbeam, ezcrab2)
    #scatter(pbeam.dist.x, pbeam.dist.px; markersize=0.1, markercolor=:red)
    RLCwake = LongitudinalRLCWake(20e9, 10.0e3, 1.0)


    oneturn=TransferMap4DChrom(opIPp, tunex, tuney, 2.0, 2.0)
    #oneturndamp=OneTurnRadiation(4000.0, 4000.0, 2000.0)
    pqbs

    line1=Lattice([ezcrab1, lb])
    line2=Lattice([ilb, ezcrab2, oneturn, lmap, RLCwake])
    longiline=Lattice([lmap, RLCwake])
end
pqbs

dpcopy=1.0.*pbeam.dist.dp
#track!(pbeam, RLCwake)

scatter(pbeam.dist.z, pbeam.dist.dp.-dpcopy)
turns=6000 #10000
coor_rec = Array{Float64, 2}(undef, 6, turns)
size_rec = Array{Float64, 2}(undef, 6, turns)
lumis = Vector{Float64}(undef, turns)
for i in 1:turns
    #track!(pbeam, line1)
    lumis[i]=track!(pbeam, estrong)
    #track!(pbeam, oneturndamp, pqbs)
    #track!(pbeam, line2)
    #track!(ebeam, longiline)
    get_emittance!(pbeam)
    coor_rec[:,i] .= pbeam.centroid
    size_rec[:,i] .= pbeam.beamsize
end

size_rec


using FFTW

plot(@view size_rec[5,:])
plot!(@view coor_rec[5,:])

fftx=fft(@view coor_rec[1,3001:end])
ffty=fft(@view coor_rec[3,3001:end])
fftz=fft(@view coor_rec[5,3001:end])
fftf=fftfreq(3000, 1.0)

#fftx=fft(@view coor_rec[1,5001:end])
#ffty=fft(@view coor_rec[3,5001:end])
#fftz=fft(@view coor_rec[5,5001:end])
#fftf=fftfreq(5000, 1.0)
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0,0.3))

get_emittance!(pbeam)
pbeam.emittance
pbeam.centroid
pbeam.beamsize


pqbs

plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=0.5)
#plot(ebeam.dist.x, ebeam.dist.px, seriestype=:scatter, markersize=0.1)
#plot(ebeam.dist.y, ebeam.dist.py, seriestype=:scatter, markersize=0.1)
#plot(ebeam., ebeam.dist.z, seriestype=:scatter, markersize=0.1)


