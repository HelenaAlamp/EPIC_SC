
struct OneTurnRadiation <: AbstractRadiation
    damping_turns::SVector{3,Float64}
    damping_coefs::SVector{3,Float64}
    damping_exp::SVector{3,Float64}
    function OneTurnRadiation(damping_turnx, damping_turny, damping_turnz)
        damping_turns=SVector{3,Float64}(damping_turnx, damping_turny, damping_turnz)
        damping_coefs=1.0./damping_turns
        damping_exp=exp.((-1.0).*damping_coefs)
        new(damping_turns, damping_coefs, damping_exp)
    end
end

function track!(ps6dcoor::AbstractVector{ps6d{T}}, or::OneTurnRadiation, eqbs::AbstractVector) where T
    ps6dcoor.x .*= or.damping_exp[1]
    ps6dcoor.px .*= or.damping_exp[1]
    ps6dcoor.y .*= or.damping_exp[2]
    ps6dcoor.py .*= or.damping_exp[2]
    ps6dcoor.z .*= or.damping_exp[3]
    ps6dcoor.dp .*= or.damping_exp[3]

    sqrt_exp=sqrt.(1 .- or.damping_exp .* or.damping_exp)
    
    ps6dcoor.x .+= (eqbs[1]*sqrt_exp[1]) .* randn.()
    ps6dcoor.px .+= (eqbs[2]*sqrt_exp[1]) .* randn.()
    ps6dcoor.y .+= (eqbs[3]*sqrt_exp[2]) .* randn.()
    ps6dcoor.py .+= (eqbs[4]*sqrt_exp[2]) .* randn.()
    ps6dcoor.z .+= (eqbs[5]*sqrt_exp[3]) .* randn.()
    ps6dcoor.dp .+= (eqbs[6]*sqrt_exp[3]) .* randn.()
    return nothing
end

function track!(beam::BunchedBeam, or::OneTurnRadiation, eqbs::AbstractVector)
    track!(beam.dist, or, eqbs)
end