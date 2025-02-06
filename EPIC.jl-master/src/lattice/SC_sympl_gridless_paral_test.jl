
function space_charge_gl_P(bb::BunchedBeam, K, Nl, Nm, a, b, Np, dt)

    gamma2lm = zeros(Nl, Nm)
    philm = zeros(Nl, Nm)

    factor = 16.0*pi*K*dt/Np/a/b


    for i in 1:Nl
        al = i * pi / a
        for j in 1:Nm
            bm = j * pi / b
            gamma2lm[i,j] = al^2 + bm^2

            temp = zeros(Threads.nthreads())
            Threads.@threads for k in 1:Np
                tid = Threads.threadid()
            
                bb.dist.x[k] += a/2
                bb.dist.y[k] += b/2
                
                temp[tid] += sin(al*bb.dist.x[k])*sin(bm*bb.dist.y[k])
            end
            
            philm[i,j] = factor*temp/gamma2lm[i,j]
            

            nthreads = Threads.nthreads()
            term1_thread = [zeros(Np) for _ in 1:nthreads]  
            term2_thread = [zeros(Np) for _ in 1:nthreads]

            
            Threads.@threads for i in 1:Nl
                tid = Threads.threadid()
                
                term1_thread[tid] -= philm[i,j]* al*cos(al * bb.dist.x[l])*sin(bm * bb.dist.y[l])
                term2_thread[tid] -= philm[i,j]* bm*sin(al * bb.dist.x[l])*cos(bm * bb.dist.y[l])

                bb.dist.px[tid] = term1_thread[tid]
                bb.dist.py[tid] = term2_thread[tid]

                bb.dist.x[l] -= a/2
                bb.dist.y[l] -= b/2


            end

        end
        
    end
    
end

