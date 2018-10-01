# Profiles transit_poly.jl
#using ForwardDiff
#using DiffResults

include("../../../src/transit_poly_struct.jl")
#include("../../../src/transit_structure.jl")

# Define function for autodiff:

#function tp_auto(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
#
## Set up call:
#  x= zeros(T,2+length(u_n))
#  x[1] = r; x[2] = b
#  for i=1:length(u_n)
#    x[2+i] = u_n[i]
#  end
#
#  # Now, define a wrapper of s2 for use with ForwardDiff:
#  function diff_tp(x::Array{T,1}) where {T <: Real}
#  # x should be a two-element vector with values [r,b]
#  rc = x[1]; bc = x[2]; uc = x[3:length(x)]
#  return transit_poly(rc,bc,uc)
#  end
#
## Now apply autodiff:
#  # Set up a type to store s_n and it's Jacobian with respect to x:
#  out = DiffResults.GradientResult(x)
#  # Compute the Jacobian (and value):
#  out = ForwardDiff.gradient!(out,diff_tp,x)
#  # Place the value in the s_2 vector:
#  tp_flux = DiffResults.value(out)
#  # And, place the Jacobian in an array:
#  tp_grad = DiffResults.gradient(out)
#return tp_flux,tp_grad
#end


#include("../../../src/compute_c_n.jl")

#include("../../../src/dJv_seriesdk.jl")
using PyPlot

# Now, run the same with transit_poly:

function profile_transit_poly(trans)
nb = 1000001; b=zeros(nb)
flux = zeros(nb); fgrad = zeros(4)
flux_grad = zeros(nb,4)
tic()
for i=0:1000000
  b[i+1] = sqrt(((i-500000.0)/500000.0*(1.0+2.0*trans.r))^2)
  trans.b = b[i+1]
  flux[i+1]= transit_poly!(trans)
  flux_grad[i+1,:] = trans.dfdrbu
end
toc()
plot(b,flux)
println("b: ",minimum(b)," ",maximum(b))
println("f: ",minimum(flux)," ",maximum(flux))
return
end

r = 0.1  # radius ratio is 0.1
u = [0.2,0.3]
b=0.0
trans = transit_init(r,b,u,true)

profile_transit_poly(trans)
