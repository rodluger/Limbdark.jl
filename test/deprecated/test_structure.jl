# Testing version which uses Transit_Struct:
using PyPlot

include("../src/transit_poly_struct.jl")

npts = 1000
iu = 30
u_n = ones(iu) / iu
r = 0.1
b0 = zeros(npts)
for i=1:npts
    b0[i] = 3.0 * (i / float(npts) - 0.5)
end

lc_ana = zeros(npts); lc_str=zeros(npts)
trans = transit_init(r,0.0,u_n,false)
for i=1:length(b0)
  b=b0[i]
  lc_ana[i] = transit_poly(r,abs(b),u_n)
  trans.b = abs(b)
  lc_str[i] = transit_poly!(trans)
  if mod(i,10) == 0
    println("i: ",i,"b: ",b," lc_ana: ",lc_ana[i]," lc_str: ",lc_str[i])
  end
end

using PyPlot

plot(b0,lc_ana, linewidth=2, label="Normal")
plot(b0,lc_str, linewidth=1, label="Struct")

