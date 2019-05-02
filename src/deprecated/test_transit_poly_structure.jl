# Tests the version of transit_poly which uses a Type
# for holding all of the variables related to the transit
# computation:

# Defines the transit Type:
#include("transit_structure.jl")
# Defines the function that computes the d_n coefficients:
#include("compute_d_n_struct.jl")

# Defines the function that computes the transit:
include("transit_poly_struct.jl")

# Instantiate a transit with r=0.1, b=0.5, u_n=[0.2,0.3]:

# Transform from u_n to d_n coefficients (this is now added to transit_init function):
#compute_d_n_grad!(trans) 

# Load in PyPlot (Julia matplotlib wrapper):
using PyPlot

# Define a function that loops over impact parameter
# and defines arrays to hold the results:

function profile_transit_poly(trans)
nb = 1000001; b=zeros(nb)
flux = zeros(nb); fgrad = zeros(4)
flux_grad = zeros(nb,4)
tic()
for i=0:1000000
  b[i+1] = sqrt(((i-500000.)/500000.*(1.+2.*trans.r))^2)
#  flux[i+1]= transit_poly!(r,b[i+1],u,fgrad)
  trans.b = b[i+1]
  flux[i+1]= transit_poly!(trans)
  flux_grad[i+1,1:2] = trans.dfdrb
  flux_grad[i+1,3:4] = trans.dfdu
end
# Now, plot the flux and derivatives:
toc()
#clf()
plot(b,flux-1,label="flux-1")
plot(b,trans.r*flux_grad[:,1],label="r*dfdr")
plot(b,trans.r*flux_grad[:,2],label="r*dfdb")
plot(b,flux_grad[:,3],label="dfdu1")
plot(b,flux_grad[:,4],label="dfdu2")
legend(loc="lower right")
println("b: ",minimum(b)," ",maximum(b))
println("f: ",minimum(flux)," ",maximum(flux))
return
end

trans = transit_init(0.1,0.5,[0.2,0.3],true)

# Call the function:
profile_transit_poly(trans)

savefig("transit_poly_derivatives.pdf", bbox_inches="tight")
