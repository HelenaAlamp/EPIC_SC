
## Beam
abstract type AbstractBeam end


## Optics
abstract type AbstractOptics end
abstract type AbstractOptics2D <:AbstractOptics end
abstract type AbstractOptics4D <:AbstractOptics end

## Elements
abstract type AbstractElement end

abstract type AbstractCavity <:AbstractElement end
abstract type AbstractAccelCavity <:AbstractCavity end
abstract type AbstractCrabCavity <:AbstractCavity end

abstract type AbstractCorrector <:AbstractElement end

abstract type AbstractDrift <:AbstractElement end

abstract type AbstractMultipole <:AbstractElement end
abstract type AbstractDipole <:AbstractMultipole end
abstract type AbstractQuadrupole <:AbstractMultipole end


abstract type AbstractWakefield <:AbstractElement end
abstract type AbstractLongiWakefield <:AbstractWakefield end
abstract type AbstractTransWakefield <:AbstractWakefield end

abstract type AbstractRadiation <: AbstractElement end

abstract type AbstractTransferMap <:AbstractElement end
abstract type AbstractTransverseMap <:AbstractTransferMap end
abstract type AbstractLongitudinalMap <:AbstractTransferMap end

abstract type AbstractLorentzBoost <:AbstractElement end

abstract type AbstractSpaceCharge <: AbstractElement end

## Lattice

abstract type AbstractLattice end
