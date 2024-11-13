mutable struct optics4DUC <: AbstractOptics4D # 4D linear transformation uncoupled
    optics_x::optics2D
    optics_y::optics2D
end
optics4DUC(bx::Float64, ax::Float64, by::Float64, ay::Float64)=optics4DUC(optics2D(bx,ax),optics2D(by,ay))

function normal_mat(op4d::optics4DUC)
    nfmx=[1.0/sqrt(op4d.optics_x.beta) 0.0; op4d.optics_x.alpha/sqrt(op4d.optics_x.beta)  sqrt(op4d.optics_x.beta)]
    nfmy=[1.0/sqrt(op4d.optics_y.beta) 0.0; op4d.optics_y.alpha/sqrt(op4d.optics_y.beta)  sqrt(op4d.optics_y.beta)] 
    return SMatrix{4,4}([nfmx zeros(2,2); zeros(2,2) nfmy])
end

function invnormal_mat(op4d::optics4DUC)
    invnfmx=[sqrt(op4d.optics_x.beta) 0.0; -op4d.optics_x.alpha*sqrt(op4d.optics_x.beta)  1.0/sqrt(op4d.optics_x.beta)]
    invnfmy=[sqrt(op4d.optics_y.beta) 0.0; -op4d.optics_y.alpha*sqrt(op4d.optics_y.beta)  1.0/sqrt(op4d.optics_y.beta)]
    return SMatrix{4,4}([invnfmx zeros(2,2); zeros(2,2) invnfmy])
end

