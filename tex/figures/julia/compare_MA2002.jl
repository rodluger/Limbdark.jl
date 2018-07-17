#include("/Users/ericagol/Courses/PlanetaryDynamics/ExoJulia/ExoJulia/Transit/occultquad.jl")
include("occultquad.jl")
include("../../../src/transit_poly.jl")

# Plot accuracy of occultquad vs. transit_poly:

nr = 10000
r=linspace(0.001,1.0,nr)
delta = 1e-12
#u1 = 0.0; u2 = 0.0
u1 = 1.0; u2 = 0.0
#u1 = 2.0; u2 = -1.0
#u1 = 0.5; u2 = 0.5

b1 = 1.0-r-delta
f_MA = occultquad.(b1,u1,u2,r)
f_AL = zeros(Float64,nr); f_AL_big = zeros(Float64,nr)
for i=1:nr
  f_AL[i] = transit_poly(r[i],b1[i],[u1;u2])
  f_AL_big[i] = convert(Float64,transit_poly(big(r[i]),big(b1[i]),[big(u1);big(u2)]))
end

b2 = r-delta
f_MA2 = occultquad.(b2,u1,u2,r);
f_AL2 = zeros(Float64,nr); f_AL_big2 = zeros(Float64,nr)
for i=1:nr
  f_AL2[i] = transit_poly(r[i],b2[i],[u1;u2])
  f_AL_big2[i] = convert(Float64,transit_poly(big(r[i]),big(b2[i]),[big(u1);big(u2)]))
end

b3 = 1.0-r+delta
f_MA3 = occultquad.(b3,u1,u2,r)
f_AL3 = zeros(Float64,nr); f_AL_big3 = zeros(Float64,nr)
for i=1:nr
  f_AL3[i] = transit_poly(r[i],b3[i],[u1;u2])
  f_AL_big3[i] = convert(Float64,transit_poly(big(r[i]),big(b3[i]),[big(u1);big(u2)]))
end

b4 = r+delta
f_MA4 = occultquad.(b4,u1,u2,r); 
f_AL4 = zeros(Float64,nr); f_AL_big4 = zeros(Float64,nr)
for i=1:nr
  f_AL4[i] = transit_poly(r[i],b4[i],[u1;u2])
  f_AL_big4[i] = convert(Float64,transit_poly(big(r[i]),big(b4[i]),[big(u1);big(u2)]))
end
using PyPlot

fig,axes = subplots(2,1)
ax = axes[1]

ax[:plot](r,f_MA,label="MA (2002), b=1-r-1e-12")
ax[:plot](r,f_AL,label="AL (2018), b=1-r-1e-12")
ax[:plot](r,f_MA2,label="MA (2002), b=r-1e-12")
ax[:plot](r,f_AL2,label="AL (2018), b=r-1e-12")
#ax[:plot](r,f_MA3,label="MA (2002), b=1-r+1e-12")
#ax[:plot](r,f_AL3,label="AL (2018), b=1-r+1e-12")
#ax[:plot](r,f_MA4,label="MA (2002), b=r+1e-12")
#ax[:plot](r,f_AL4,label="AL (2018), b=r+1e-12")
ax[:legend](loc="lower left",fontsize=8)
ax[:set_xlabel]("r; b = 1-r, b=r")
ax[:set_ylabel]("Flux")
ax[:axis]([-0.01,1.01,0,1.01])

ax = axes[2]
ax[:semilogy](r,abs.(f_MA-f_AL_big),label="MA (2002), b=1-r-1e-12")
ax[:semilogy](r,abs.(f_AL-f_AL_big),label="AL (2018), b=1-r-1e-12")
ax[:semilogy](r,abs.(f_MA2-f_AL_big2),label="MA (2002), b=r-1e-12")
ax[:semilogy](r,abs.(f_AL2-f_AL_big2),label="AL (2018), b=r-1e-12")
#ax[:semilogy](r,abs.(f_MA3-f_AL_big3),label="MA (2002), b=1-r+1e-12")
#ax[:semilogy](r,abs.(f_AL3-f_AL_big3),label="AL (2018), b=1-r+1e-12")
#ax[:semilogy](r,abs.(f_MA4-f_AL_big4),label="MA (2002), b=r+1e-12")
#ax[:semilogy](r,abs.(f_AL4-f_AL_big4),label="AL (2018), b=r+1e-12")
ax[:set_xlabel]("r; b= 1-r, b=r")
ax[:set_ylabel]("Error")
#ax[:legend](loc="left",fontsize=8)
ax[:axis]([-0.01,1.01,1e-18,1e-4])

savefig("compare_MA2002_linear.pdf", bbox_inches="tight")
