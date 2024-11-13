
using Revise, EPIC
using Test
using BenchmarkTools
using Plots

theme(:ggplot2)

begin
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    pbeam=BunchedBeam(PROTON, 0.688e11, 275e9,  1000, [11.3e-9, 1e-9, 3.7e-2])
    #initilize_zsliceSC!(pbeamgauss, :gaussian, :evennpar, 5.0)
    mainRF=AccelCavity(591e6, 15.8e6, 7560.0, π)
    αc=1.5e-3
    lmap=LongitudinalRFMap(αc, mainRF)   
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)
end

#useless for SC
get_2nd_moment!(pbeam)
get_centroid!(pbeam)
get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize

SC_point!(pbeam)

#fieldvecSC = [0.0 , 0.0, 0.0]
#Bassetti_ErskineSC!(fieldvecSC, 0.0, 0.0, pbeam.beamsize[1], pbeam.beamsize[3])

#plot:

plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=0.5)
plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=0.5)


#repeat for number of turns (only SC)
#= 
turns = 1000
for i in 3:turns
    if i%3 == 0
    get_emittance!(pbeam)
    coor_rec[:,i] .= pbeam.centroid
    size_rec[:,i] .= pbeam.beamsize
    end
end=#

