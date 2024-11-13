module EPIC
using StaticArrays
using StructArrays

include("abstypes.jl")
include("beam/beam.jl")

include("lattice/elements.jl")
include("lattice/lattice.jl")
include("optics/optics.jl")

export PROTON, ELECTRON, GOLDION
export ps6d

export BunchedBeam
export get_centroid!,get_2nd_moment!, get_emittance!
export initilize_6DGaussiandist!, histogram1DinZ!

export optics4DUC

export StrongGaussianBeam, StrongThinGaussianBeam, initilize_zslice!, Bassetti_Erskine!, crab_crossing_setup!, get_tune_spread!

export SC_lens, SC_kick!

export Drift, ThinCorrector
export TransferMap4D, Inverse_TransferMap4D, TransferMap4DChrom, LongitudinalRFMap, OneTurnRadiation
export AccelCavity, CrabCavity, easyCrabCavity, LongitudinalRLCWake
export LorentzBoost, InvLorentzBoost
export Lattice
export track!
export get_synchrotron_tune

#export shape_function, d_shape_function, calculate_K
export calculate_K, space_charge_gl, space_charge_gl_P
export SPACECHARGE, SC_gl_track!


end
