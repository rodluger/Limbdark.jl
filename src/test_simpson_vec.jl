
include("deprecated/simpson_vec.jl")

function f(x,f_of_x)
f_of_x = [x,sin(x),cos(x),x^3]
return
end

#f(x) = 

i = 0.0

I_of_f = zeros(3)

eps = 1e-8

simpson_vec(0.0,float(pi),f,I_of_f,i,eps,8,3) #-[2.0,0.0,pi^4/4]

