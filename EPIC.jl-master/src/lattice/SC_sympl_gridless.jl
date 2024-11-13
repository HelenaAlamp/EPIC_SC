
#calculate generalized perveance / current #total charge for PROTON: 1.1023e-8
function calculate_K(bb::BunchedBeam, σz)

    K = 2*bb.particle.classrad0*bb.num_particle/(bb.beta^3*bb.gamma^3*σz)
    
    return K
end


function space_charge_gl(bb::BunchedBeam, K, Nl, Nm, a, b, Np, dt)

    gamma2lm = zeros(Nl, Nm)
    philm = zeros(Nl, Nm)

    factor = 16.0*pi*K*dt/Np/a/b

    #x = zeros(Np)
    #y = zeros(Np)

    for i in 1:Nl
        al = i * pi / a
        for j in 1:Nm
            bm = j * pi / b
            gamma2lm[i,j] = al^2 + bm^2

            temp = 0.0
            for k in 1:Np
                bb.dist.x[k] += a/2
                bb.dist.y[k] += b/2
                temp += sin(al*bb.dist.x[k])*sin(bm*bb.dist.y[k])

            end
            philm[i,j] = factor*temp/gamma2lm[i,j]
            
            for l in 1:Np
                bb.dist.px[l] -= philm[i,j]* al*cos(al * bb.dist.x[l])*sin(bm * bb.dist.y[l])
                bb.dist.py[l] -= philm[i,j]* bm*sin(al * bb.dist.x[l])*cos(bm * bb.dist.y[l])

                bb.dist.x[l] -= a/2
                bb.dist.y[l] -= b/2

            end

        end
        
    end
    
end


mutable struct SPACECHARGE <: AbstractSpaceCharge
    effective_len::Float64
    Nl::Int64
    Nm::Int64
    a::Float64
    b::Float64
    eletype::String

    function SPACECHARGE(effective_len::Float64, Nl::Int64, 
                        Nm::Int64, a::Float64, b::Float64)
        new(effective_len, Nl, Nm, a, b) #Nl,m = a,b/num_particle
    end
end


function SC_gl_track!(ele::SPACECHARGE, bb::BunchedBeam, num_particles::Int64)
   
    get_emittance!(bb)
    σz = bb.beamsize[5]

    K = calculate_K(bb, σz)

    ele.a = bb.beamsize[1]*10
    ele.b = bb.beamsize[3]*10
    
    space_charge_gl(bb, K, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, ele.effective_len)
    
end


########################### test ###########################

###### test calculate_K ######
#=
num_particles=1000
beam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-3])
K = calculate_K(beam,0.006)
beam.dist.x

###### for test ######
#generate beam distribution
opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0) 
mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
αc=1.5e-3
lmap = LongitudinalRFMap(αc, mainRF)
initilize_6DGaussiandist!(beam, opIPp, lmap) #vector 6D phase space


Nl = 15
Nm = 15
#dx = 3e-5
#dy = 3e-5
a = 10e-3 #beam pipe 10mm
b = 10e-3
Np = num_particles
dt = 0.5 #3800/5/3e8 #ring division - time length

sc = SPACECHARGE(dt, Nl, Nm, a, b)

space_charge_gl(beam, K, Nl, Nm, a, b, Np, dt)

SC_gl_track!(sc, beam, Np)

using Plots
plot(beam.dist.x, px.- px2, seriestype=:scatter, markersize=2)
=#
