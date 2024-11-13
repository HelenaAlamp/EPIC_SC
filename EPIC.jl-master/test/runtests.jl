using Revise, EPIC
using Test
using BenchmarkTools
using StructArrays
using Plots
using LinearAlgebra
using StaticArrays
using SpecialFunctions

using Distributions

# test bunched beam with IP optics and 6D Gaussian distribution

pbeam=BunchedBeam(PROTON, 1e11, 250e9,  1000000, [2.5e-6, 1e-6, 1e-1])

opIP=optics4DUC(4.0,0.0,1.0,0.0)

mainRF=AccelCavity(197e6, 1e6, 2520.0, 0.0)
initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3)
minimum(pbeam.dist.y), maximum(pbeam.dist.y)
@benchmark begin
    initilize_6DGaussiandist!(pbeam, opIP, mainRF, 1.8e-3)
    #track!(pbeam.dist, mainRF, 1.8e-3)
    #get_centroid!(pbeam)
    
end

get_2nd_moment!(pbeam)

get_centroid!(pbeam)
@btime begin get_centroid!(pbeam) end
get_emittance!(pbeam)
pbeam.emittance
sqrt(pbeam.moment2nd[5,5])

# test drift, checked efficiency
drift1=Drift(1.0)
track!(pbeam, drift1) 
@btime begin track!(pbeam, drift1) end

# test corrector, checked efficiency
corrector1=ThinCorrector(1e-3, 1e-3)
track!(pbeam.dist, corrector1)
@btime begin track!(pbeam.dist, corrector1) end

# test accRF, checked efficiency
track!(pbeam, mainRF)
@btime begin 
    track!(pbeam, mainRF) 
    track!(pbeam, mainRF)
    track!(pbeam, mainRF)
end

# test crab CrabCavity
crab1=CrabCavity(400e6, 1e6)
track!(pbeam, crab1)
@btime begin track!(pbeam, crab1) end

# test Strong Beam setup
sigz=0.1
pstbeam=StrongGaussianBeam(PROTON, 1e11, 250e9,  opIP, [4e-6, 1e-6, sigz], 20)
initilize_zslice!(pstbeam, :gaussian, :evennpar, 100.0)
pstbeam.zslice_center

# test histogram and wakefield
pbeam.znbin

histogram1DinZ!(pbeam)
minimum(pbeam.inzindex)
pbeam.zhist_edges
pbeam.zhist
zsep=(pbeam.zhist_edges[end]-pbeam.zhist_edges[1])/pbeam.znbin
@benchmark begin
    histogram1DinZ!(pbeam)
end

#zhist, zhist_edges=histogram1DinZ!(pbeam, 200)
plot(pbeam.zhist_edges[1:end-1], pbeam.zhist, seriestype=:scatter, markersize=1, legend=false, xlabel="z [m]", ylabel="N", title="6D Gaussian distribution")




RLCwake = LongitudinalRLCWake(1e9, 1e5, 1.0)
old_dp=pbeam.dist.dp.*1.0
pbeam.ztemp4
track!(pbeam, RLCwake)


@btime begin track!(pbeam, RLCwake) end

plot(pbeam.dist.z, pbeam.dist.dp.-old_dp, markersize=1, seriestype=:scatter, legend=false, xlabel="z [m]", ylabel="dp/p" )


## test transfer map
oneturn=TransferMap4DChrom(opIP, 0.08, 0.139, 1.0, 1.0)
oneturn.invumat
track!(pbeam, oneturn)
@benchmark begin track!(pbeam, oneturn) end


oneturnrad=OneTurnRadiation(4000.0, 4000.0, 2000.0)
track!(pbeam, oneturnrad, [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5])
@benchmark begin track!(pbeam, oneturnrad, [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]) end

# test B-E formula
fieldvec = MVector{3, Float64}(undef)
fieldvec[1]

@benchmark begin pbeam end

@benchmark begin for i in 1:1 Bassetti_Erskine!(fieldvec, pbeam.dist.x[1], pbeam.dist.y[1], 1e-4, 1e-5) end end
fieldvec
Bassetti_Erskine(fieldvec, -3e-4, 3e-5, 1e-4, 1e-5)
fieldvec




# test Weak-Strong Simulation

pbeam=BunchedBeam(PROTON, 0.688e11, 275e9,  1000000, [11.3e-9, 1e-9, 3.75e-2])
opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
mainRF=AccelCavity(591e6, 19.6e6, 7560.0, 0.0)
initilize_6DGaussiandist!(pbeam, opIPp, mainRF, 1.8e-3)
get_emittance!(pbeam)
pbeam.beamsize

opIPe=optics4DUC(0.45, 0.0, 0.056, 0.0)
estrong=StrongGaussianBeam(ELECTRON, 1.72e11, 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 1)

initilize_zslice!(estrong, :gaussian, :evennpar, 10.0)
#lumi=track!(pbeam, estrong) 
#lumi*3e8/3833*1320*0.9/1e4
#pbeam.dist.x

#estrong.zslice_center
#estrong.zslice_npar

@benchmark begin
    track!(pbeam, estrong)
end
a=randn(10)
b=@view a[1:end-1] .+ (@view a[2:end])



line=Lattice()
pbeam=BunchedBeam(PROTON, 0.688e11, 275e9,  1000000, [11.3e-9, 1e-9, 3.75e-2])
ebeam=BunchedBeam(ELECTRON, 1.72e11, 10e9,  100000, [20e-9, 1.3e-9, 3.75e-3])

track!(pbeam, line)

@testset "EPIC.jl" begin
    # Write your tests here.
    a=1
end


