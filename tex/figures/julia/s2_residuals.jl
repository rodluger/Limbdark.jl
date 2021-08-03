# Uses new formulation from limbdark paper.
using Limbdark
import Limbdark: Transit_Struct, s2

#include("s2_stable.jl")
#include("../../../src/s2.jl")
#include("../../../src/define_constants.jl")

# Run some tests:

using PyPlot
nb = 1001; nr = 1001
b=range(0.0,stop=2,length=nr)
r=range(0.0,stop=2,length=nr)
fgrid = zeros(Float64,nr,nb)
fgrid_old = zeros(Float64,nr,nb)
fgrid_big = zeros(Float64,nr,nb)
for ib=1:nb
#  if mod(ib,10) == 0
#    println("Finished: ",ib/nr*100,"%")
#  end
  for ir=1:nr
    fgrid[ir,ib]=s2(r[ir],b[ib])
    fgrid_big[ir,ib]=convert(Float64,s2(big(r[ir]),big(b[ib])))
  end
end

clf()
figure(figsize=(5, 4))
img=imshow(1e15 * (fgrid-fgrid_big),interpolation="nearest",origin="lower",
           extent=[0,2,0,2],cmap="plasma",vmin=-1.5, vmax=1.5)
plot([0,1],[1,2],"k-")
plot([1,2],[0,1],"k-")
plot([0,0.38],[0,0.38],"k-")
plot([0.5,2],[0.5,2],"k-")
plot([0,1],[1,0],"k-")
annotate(L"$r = 1 + b$", xy=(0.5, 1.5), xytext=(0, 0), xycoords="data",
         textcoords="offset points", ha="center", va="bottom",
         fontsize=12, color="k", rotation=45)

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
ticks = [-1.5, -1, -0.5, 0, 0.5, 1, 1.5]
cbar = colorbar(img, ticks=ticks)
cbar[:ax][:set_yticklabels]([L"$-1.5\times 10^{-15}$",
                             L"$-1.0\times 10^{-15}$",
                             L"$-0.5\times 10^{-15}$",
                             L"$0$",
                             L"$0.5\times 10^{-15}$",
                             L"$1.0\times 10^{-15}$",
                             L"$1.5\times 10^{-15}$"])

savefig("s2_residuals.pdf")
