

# Computing light curves errors as a function of the order of the limb-darkening
include("../../../test/loglinspace.jl")
#include("../src/transit_poly_struct.jl")
#include("../src/nonlinear_limbdarkening/get_coeff.jl")

import Limbdark

using PyPlot
using Statistics

n_u_max =  30
r = 0.1
b = 0.01
nb = 1000
b0 = 0.01
bgrid = sqrt.(b0^2 .+linearspace(-1.2,1.2,nb).^2)
fgrid = zeros(nb)
fbig  = zeros(BigFloat,nb)
ntrial = 10
err_vs_n = zeros(ntrial,n_u_max)
# Use the "power-two" law to get the coefficients:
#nmu = 100; mu = logspace(-big(3.0),big(0.0),nmu)
#alpha = 1.25
fig,axes = subplots(1,2)
for j=1:ntrial
for n_u=1:n_u_max
  #u = ones(n_u)/n_u
  u = rand(n_u)
  u ./= sum(u)
  t = Limbdark.transit_init(r,b,u,false)
  # Pick the coefficients for the non-linear limb-darkening law:
#  fnon = mu.^alpha-mu
#  ci,fmod = get_coeff!(mu,fnon,n_u,nmu)
#  ax = axes[3]
#  ax[:plot](mu,fnon)
#  ax[:plot](mu,fmod,".")
#  t.g_n[1] = 0.0
#  t.g_n[2] = 1-2sum(ci)
#  t.g_n[3:n_u+1]=ci
#  t.den = 1/(pi*(t.g_n[1]+2//3*t.g_n[2]))
  t_big = Limbdark.transit_init(big(r),big(b),big.(u),false) # BigFloat
#  t_big.g_n[1] = big(0.0)
#  t_big.g_n[2] = 1-2sum(ci)
#  t_big.g_n[3:n_u+1]=ci
#  t_big.den = 1/(pi*(t_big.g_n[1]+2//3*t_big.g_n[2]))
  # See if higher accuracy g_n values helps:
  # t.g_n = convert(Array{Float64,1},t_big.g_n)
  for i=1:nb
    t.b = bgrid[i]
    fgrid[i] = Limbdark.transit_poly_g(t)
    t_big.b = big(bgrid[i])
    fbig[i] = Limbdark.transit_poly_g(t_big)
  end
  logerr = (fgrid .-convert(Array{Float64,1},fbig))/maximum(1.0 .-fgrid)
  ax = axes[1]
  ax[:semilogy](bgrid,logerr)
  ax = axes[2]
  ax[:semilogy](bgrid,fgrid)
  println("n_u: ",n_u," max error: ",maximum(logerr))
  err_vs_n[j,n_u] = maximum(logerr)
#  read(STDIN,Char)
end
end
#read(STDIN,Char)
clf()
err_med = zeros(n_u_max)
for j=1:ntrial
  semilogy(err_vs_n[j,:],".")
end
for i=1:n_u_max
  err_med[i] = median(err_vs_n[:,i])
end
#semilogy(err_vs_n)
semilogy(err_med,linewidth=3)
xlabel("Order of limb-darkening")
ylabel("Fractional error relative to maximum transit depth")
title(L"Errors on transit computation versus order of limb-darkening for $r=0.1$")
savefig("fractional_error_vs_order.pdf", bbox_inches="tight")

