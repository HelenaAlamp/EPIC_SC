struct Drift <: AbstractDrift
    length::Float64
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, et, normlen, oo::Drift) where T
    et .= (one(T) .+ ps6dcoor.dp)
    normlen .= oo.length ./ (et .* et .- ps6dcoor.px .* ps6dcoor.px .- ps6dcoor.py .* ps6dcoor.py)
    ps6dcoor.x .+= normlen .* ps6dcoor.px 
    ps6dcoor.y .+= normlen .* ps6dcoor.py 
    ps6dcoor.z .+=  et .* normlen .- oo.length
    nothing
end

function track!(beam::BunchedBeam, oo::Drift) 
    track!(beam.dist, beam.temp1, beam.temp2, oo)
end