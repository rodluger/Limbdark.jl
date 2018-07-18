# Profiles transit_poly.jl

include("../../../src/transit_poly_struct.jl")
#include("../../../src/transit_structure.jl")

r = 0.1  # radius ratio is 0.1
u = [0.2,0.3]
b=0.0
trans = transit_init(r,b,u,true)

#include("../../../src/compute_c_n.jl")

#include("../../../src/dJv_seriesdk.jl")
using PyPlot

# Now, run the same with transit_poly:

#function profile_transit_poly()
nb = 1000001; b=zeros(nb)
flux = zeros(nb); fgrad = zeros(4)
flux_big = zeros(nb); fgrad_big = zeros(BigFloat,4)
flux_grad = zeros(nb,4)
flux_grad_big = zeros(nb,4)
tic()
for i=0:1000000
  b[i+1] = sqrt(((i-500000.)/500000.*(1.+2.*r))^2)
  trans.b = b[i+1]
  flux[i+1]= transit_poly!(trans)
  flux_grad[i+1,:] = trans.dfdrbu
end
toc()
plot(b,flux)
println("b: ",minimum(b)," ",maximum(b))
println("f: ",minimum(flux)," ",maximum(flux))
#return
#end
