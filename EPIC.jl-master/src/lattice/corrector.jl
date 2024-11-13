struct ThinCorrector <: AbstractCorrector
    dx::Float64
    dy::Float64
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, tc::ThinCorrector) where T
    ps6dcoor.px .+= tc.dx
    ps6dcoor.py .+= tc.dy
    nothing
end

function track!(beam::BunchedBeam, tc::ThinCorrector)
    track!(beam.dist, tc)
    nothing
end