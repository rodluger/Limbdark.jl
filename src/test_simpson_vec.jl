include("simpson_vec.jl")

f(x) = [sin(x),cos(x),x^3]

i = 0.0

I_of_f = zeros(3)

eps = 1e-8

simpson_vec(0.0,float(pi),f,I_of_f,i,eps,8,3)-[2.0,0.0,pi^4/4]

