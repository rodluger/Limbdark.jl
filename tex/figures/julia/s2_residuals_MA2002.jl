# Uses new formulation from limbdark paper.
using Limbdark
import Limbdark: Transit_Struct, s2

#include("s2_stable.jl")
#include("../../../src/s2.jl")
#include("../../../src/define_constants.jl")

if VERSION >= v"0.7"
  using DelimitedFiles
end

using PyPlot
nb = 1001; nr = 1001
b=range(0.0,stop=2,length=nr)
r=range(0.0,stop=2,length=nr)

# We compute the `occultquad` flux in
# a separate Python call
fgrid = readdlm("s2_MA2002.txt")

fgrid_old = zeros(Float64,nr,nb)
fgrid_big = zeros(Float64,nr,nb)
resid = zeros(Float64,nr,nb)
for ib=1:nb
  for ir=1:nr
    fgrid_big[ir,ib]=convert(Float64,s2(big(r[ir]),big(b[ib])))
    resid[ir,ib] = abs.(fgrid[ir,ib]-fgrid_big[ir,ib])
    if resid[ir,ib] < 1e-15
      resid[ir,ib] = 1e-15
    end
  end
end

clf()
figure(figsize=(5, 4))
logresid = log10.(resid)
img=imshow(logresid,interpolation="nearest",origin="lower",
           extent=[0,2,0,2],cmap="plasma",vmin=-15, vmax=-3)
plot([0,1],[1,2],"k-")
plot([1,2],[0,1],"k-")
plot([0,0.38],[0,0.38],"k-")
plot([0.5,2],[0.5,2],"k-")
plot([0,1],[1,0],"k-")
annotate(L"$r = 1 + b$", xy=(0.5, 1.5), xytext=(0, 0), xycoords="data",
         textcoords="offset points", ha="center", va="bottom",
         fontsize=12, color="w", rotation=45)

annotate(L"$r = b$", xy=(1, 1), xytext=(0, 0), xycoords="data",
         textcoords="offset points", ha="center", va="bottom",
         fontsize=12, color="k", rotation=45)

annotate(L"$r = b - 1$", xy=(1.5, 0.5), xytext=(0, 0), xycoords="data",
         textcoords="offset points", ha="center", va="bottom",
         fontsize=12, color="k", rotation=45)

annotate(L"$r = 1 - b$", xy=(0.5, 0.5), xytext=(-9, -9), xycoords="data",
         textcoords="offset points", ha="center", va="center",
         fontsize=12, color="k", rotation=-45)

xlabel("Impact parameter")
ylabel("Radius ratio")
ticks = [-15, -12, -9, -6, -3]
cbar = colorbar(img, ticks=ticks)
cbar[:ax][:set_yticklabels]([L"$10^{-15}$",
                             L"$10^{-12}$",
                             L"$10^{-9}$",
                             L"$10^{-6}$",
                             L"$10^{-3}$"])

savefig("s2_residuals_MA2002.pdf")
