#test for inverse transfer matrix - for transfermap.jl
#=
using StaticArrays
using LinearAlgebra
#using Base

tunex=0.228
tuney=0.21
s1,c1=sincos(2π*tunex)
r1=@SMatrix [c1 s1; -s1 c1]
s2,c2=sincos(2π*tuney)
r2=@SMatrix [c2 s2; -s2 c2]
zm=@SMatrix [0.0 0.0; 0.0 0.0]
rotation=SMatrix{4,4}([r1 zm; zm r2])
tune=SVector{2,Float64}(tunex, tuney)

inv_rotation = inv(rotation)

println(rotation)
println(inv_rotation)

#devo provare con la matrice S intera (invnormal_mat(o2)*rotation*normal_mat(o1))

#=
s1,c1=sincos(2π*(tunex))
r1=@SMatrix [c1 -s1; s1 c1]
s2,c2=sincos(2π*(tuney))
r2=@SMatrix [c2 -s2; s2 c2]
zm=@SMatrix [0.0 0.0; 0.0 0.0]
rotation=SMatrix{4,4}([r1 zm; zm r2])
tune=SVector{2,Float64}(tunex, tuney)

println(rotation)

if rotation == inv_rotation
    println("OK")
else println("no")
end
=#

# check inv(M)*M = I

M = inv_rotation*rotation
M = MArray(M)
threshold = 1e-16   #parameter threshold

# condition to remove calculation error for inv(M)
for i in eachindex(M)
    if abs(M[i]) < threshold
        M[i] = 0
    end
end
 

#function to remove calculation error for inv(M)
function threshold_to_zero!(matrix::Matrix, threshold::Float64)
    # Vectorize the threshold comparison and assignment
    matrix .= (matrix > threshold) .* matrix
end

#apply function
threshold_to_zero!(M, 1e-10)


#test on det theorem det(inv(M)) = 1/det(M)
if det(inv_rotation) == 1/det(rotation)
    println("OK")
else println("no")
end
=#