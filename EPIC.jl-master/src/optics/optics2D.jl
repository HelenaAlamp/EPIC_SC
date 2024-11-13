

mutable struct optics2D <: AbstractOptics2D
    beta::Float64
    alpha::Float64
    gamma::Float64
    phase::Float64
    eta::Float64
    etap::Float64
end
optics2D(b::Float64, a::Float64, ph::Float64, eta::Float64, etap::Float64)=optics2D(b,a,(1+a^2)/b,ph,eta,etap)
optics2D(b::Float64, a::Float64)=optics2D(b,a,(1+a^2)/b,0.0,0.0,0.0)

function normal_mat(op2d::optics2D)
    return @SMatrix [1.0/sqrt(op2d.beta) 0.0; op2d.alpha/sqrt(op2d.beta)  sqrt(op2d.beta)]
end

function invnormal_mat(op2d::optics2D)
    return @SMatrix [sqrt(op2d.beta) 0.0; -op2d.alpha*sqrt(op2d.beta)  1.0/sqrt(op2d.beta)]
end


