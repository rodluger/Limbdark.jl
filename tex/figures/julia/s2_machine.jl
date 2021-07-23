# Uses new formulation from limbdark paper
using Limbdark
import Limbdark: Transit_Struct, s2

#include("s2_stable.jl")
#include("../../../src/transit_structure.jl")
include("../../../test/loglinspace.jl")
#include("../../../src/s2.jl")
#include("../../../src/define_constants.jl")

s2diff(b,r) = s2(r,b)-convert(Float64,s2(big(r),big(b)))

# Run some tests:

using PyPlot

# Found solution for using mathfrak here:
# https://github.com/JuliaPlots/Plots.jl/issues/997
#
using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["text.latex.preamble"] = ["\\usepackage{eufrak}",]

nb = 1024; nr = 1024
fig,axes = subplots(3,2, figsize=(13, 8))

# Compute r = b+-\epsilon
epsilon = 1e-8
fepsp= zeros(Float64,nb)
fepsbigp= zeros(Float64,nb)
fepsm= zeros(Float64,nb)
fepsbigm= zeros(Float64,nb)
#b=range(epsilon,stop=2.0,length=nr)
b=linearspace(epsilon,2.0,nr)
for ib=1:nb
  fepsbigp[ib] = float(s2(big(b[ib]+epsilon),big(b[ib])))
  fepsp[ib] = float(s2(b[ib]+epsilon,b[ib]))
  fepsbigm[ib] = float(s2(big(b[ib]-epsilon),big(b[ib])))
  fepsm[ib] = float(s2(b[ib]-epsilon,b[ib]))
end

ax = axes[1]
ax[:plot](b,s2.(b,b),"-", lw=2,color="k",alpha=0.3,label=L"$r=b$")
ax[:plot](b,fepsp,"--", lw=2,color="C0",label=L"$r=b+\epsilon$")
ax[:plot](b,fepsm,":", lw=2,color="C1",label=L"$r=b-\epsilon$")
ax[:set_ylabel](L"$\mathfrak{s}_1$", fontsize=14)
ax[:set_title](L"$\mathfrak{s}_1(r = b \pm \epsilon,b)$", fontsize=16)
ax[:legend](loc="upper right",fontsize=10)

ax = axes[2]
ax[:plot](b,(fepsp-s2.(b,b))/epsilon,lw=2,color="C0",label=L"$\mathrm{double}$")
ax[:plot](b,(fepsbigp-s2.(b,b))/epsilon,linewidth=2,linestyle="--", color="C1",label=L"$\mathrm{double}$")
ax[:plot](b,(fepsm-s2.(b,b))/epsilon,lw=2,color="C0")
ax[:plot](b,(fepsbigm-s2.(b,b))/epsilon,linewidth=2,linestyle="--",color="C1")
ax[:set_ylabel](L"$\Delta \mathfrak{s}_1 / \epsilon$", fontsize=14)

ax[:annotate](L"$r = b - \epsilon$", xy=(1.8, 1.8), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=0)

ax[:annotate](L"$r = b + \epsilon$", xy=(1.8, -1.8), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=0)

ax = axes[3]
ax[:plot](b,1e15 * s2diff.(b,b .+ epsilon),label=L"$r=b+\epsilon$")
ax[:plot](b,1e15 * s2diff.(b,b .- epsilon),linewidth=2,linestyle="--",label=L"$r=b-\epsilon$")
ax[:set_ylabel](L"$\mathrm{error}\ \times 10^{15}$", fontsize=14)
ax[:set_xlabel](L"$b$", fontsize=14)
ax[:legend](loc="lower right",fontsize=10)

# Compute b+r = 1+-\epsilon
epsilon = 1e-8
#b=range(0.0,stop=2.0,length=nr)
b=linearspace(0.0,2.0,nr)
for ib=1:nr
  fepsbigp[ib] = float(s2(big(1.0 .- b[ib] .+ epsilon),big(b[ib])))
  fepsp[ib] = float(s2(1.0 .- b[ib] .+ epsilon,b[ib]))
  fepsbigm[ib] = float(s2(big(1.0 .- b[ib] .- epsilon),big(b[ib])))
  fepsm[ib] = float(s2(1.0 .- b[ib] .- epsilon,b[ib]))
end

ax = axes[4]
ax[:plot](b,s2.(1.0 .- b,b),"k-",label=L"$r=1-b$",lw=2,alpha=0.3)
ax[:plot](b,fepsp,"--",label=L"$r=1-b+\epsilon$",lw=2,color="C0")
ax[:plot](b,fepsm,":",label=L"$r=1-b-\epsilon$",lw=2,color="C1")
ax[:legend](loc="lower right",fontsize=10)
ax[:set_title](L"$\mathfrak{s}_1(r = 1 - b \pm \epsilon,b)$", fontsize=16)

ax = axes[5]
ax[:plot](b,(fepsp-s2.(1.0 .- b,b))/epsilon,lw=2,label=L"$\mathrm{double}$", color="C0")
ax[:plot](b,(fepsbigp-s2.(1.0 .- b,b))/epsilon,lw=2,linestyle="--",label=L"$\mathrm{BigFloat}$", color="C1")
ax[:plot](b,(fepsm-s2.(1.0 .- b,b))/epsilon,lw=2, color="C0")
ax[:plot](b,(fepsbigm-s2.(1.0 .- b,b))/epsilon,lw=2,linestyle="--", color="C1")
ax[:legend](loc="lower right",fontsize=10)

ax[:annotate](L"$r = 1-b - \epsilon$", xy=(0.75, 1.5), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=-35)

ax[:annotate](L"$r = 1-b + \epsilon$", xy=(0.75, -1.5), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation= 35)

ax = axes[6]
ax[:plot](b,1e15*s2diff.(1.0 .- b .+ epsilon,b),label=L"$r=1-b+\epsilon$")
ax[:plot](b,1e15*s2diff.(1.0 .- b .- epsilon,b),linewidth=2,linestyle="--",label=L"$r=1-b-\epsilon$")
ax[:legend](loc="lower right",fontsize=10)
ax[:set_xlabel](L"$b$", fontsize=14)

savefig("s2_machine.pdf", bbox_inches="tight")
