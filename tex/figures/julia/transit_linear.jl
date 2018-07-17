# Uses new formulation from limbdark paper.

include("../../../src/s2.jl")

# Run some tests:

using PyPlot
nb = 401; nr = 401
b=linspace(0.0,2,nr)
r=linspace(0.0,2,nr)
fgrid = zeros(Float64,nr,nb)
for ib=1:nb
#  if mod(ib,10) == 0
#    println("Finished: ",ib/nr*100,"%")
#  end
  for ir=1:nr
    fgrid[ir,ib],Eofk,Em1mKdm=s2(r[ir],b[ib])
  end
end

clf()
img=imshow(fgrid,interpolation="nearest",origin="lower",extent=[0,2,0,2])
plot([0,1],[1,2],label=L"$r=1+b$")
plot([1,2],[0,1],label=L"$r=b-1$")
plot([0,2],[0,2],label=L"$r=b$")
plot([0,1],[1,0],label=L"$r=1-b$")
legend(loc="upper left",fontsize=10)
xlabel("Impact parameter")
ylabel("Radius ratio")
colorbar(img)
savefig("transit_linear.pdf", bbox_inches="tight")
