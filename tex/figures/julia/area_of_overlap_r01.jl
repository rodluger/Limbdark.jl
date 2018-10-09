function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}
# How to compute (twice) area squared of triangle with
# high precision (Goldberg 1991):
a,b,c=reverse(sort([a,b,c]))
area = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
return area
end

# Shows the overlap area accuracy for the two formulae:

# 1). Mandel & Agol:

function circle_overlap(r,b)
kap1 = acos((1-r^2+b^2)/(2*b))
kap0 = acos((r^2+b^2-1)/(2*b*r))
akite = sqrt((4*b^2-(1+b^2-r^2)^2)/4)
area = r^2*kap0+kap1-akite

# 2). Luger & Agol:

kite_area2 = sqrt(sqarea_triangle(one(r),b,r))
# Angle of section for occultor:
kap0  = atan2(kite_area2,(r-1)*(r+1)+b^2)
# Angle of section for source:
kap1 = atan2(kite_area2,-((r-1)*(r+1)-b^2))
# Flux of visible uniform disk:
sn = kap1 + r^2*kap0 - .5*kite_area2

# Now at high precision:
r_big = big(r); b_big=big(b)
kite_area2_big = sqrt(sqarea_triangle(big(1.0),b_big,r_big))
# Angle of section for occultor:
kap0_big  = atan2(kite_area2,(r_big-1)*(r_big+1)+b_big^2)
# Angle of section for source:
kap1_big = atan2(kite_area2_big,-((r_big-1)*(r_big+1)-b_big^2))
# Flux of visible uniform disk:
sn_big = kap1_big + r_big^2*kap0_big - 0.5*kite_area2_big

return area,sn,sn_big
end


using PyPlot
fig,axes = subplots(2,2)

r = 0.1
nb = 1000
db = logspace(-15,-2,nb)
b = 1-r+db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(BigFloat,nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[1]
ax[:plot](db,abs.(pi*r^2-area),label="r=0.1, MA(2002)",lw=3)
ax[:plot](db,abs.(pi*r^2-sn),label="r=0.1, AL(2018)",".")
ax[:plot](db,abs.(convert(Array{Float64,1},pi*big(r)^2-sn_big)),label="r=0.1, BigFloat",linestyle="--",lw=3)
#ax[:set_xlabel](L"$\log_{10}(b-(1-r))$")
ax[:set_title ](L"$\pi r^2$-Area of overlap")
ax[:legend](loc="lower right",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-22,1e-2])
t1 = linspace(0,2pi,100)
aspect = 1.4*24.0/17.0
x1 = -12.0+cos.(t1); y1=-5.0+aspect*sin.(t1)
x2 = -11.1+0.1*cos.(t1); y2=-5.0+0.1*aspect*sin.(t1)
ax[:plot](x1,y1)
ax[:plot](x2,y2)
ax[:plot]([1e-16,1e-2],[1,1]*(pi*0.1^2/2^53),c="grey",linestyle="-.")

ax = axes[2]
ax[:plot](db,abs.(area-convert(Array{Float64,1},sn_big)),label="r=0.1, MA(2002)",lw=3)
ax[:plot](db,abs.(sn-convert(Array{Float64,1},sn_big)),label="r=0.1, AL(2018)",".")
ax[:set_xlabel](L"$b-(1-r)$")
ax[:set_ylabel]("Abs(Error in area of overlap)")
ax[:legend](loc="center left",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-25,1e-5])
ax[:plot]([1e-16,1e-2],[1,1]*(pi*0.1^2/2^53),c="grey",linestyle="-.")

b = 1+r-db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(BigFloat,nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[3]
mask = area .> 0.0
ax[:plot](db[mask],area[mask],label="r=0.1, MA(2002)",lw=3)
ax[:plot](db,sn,label="r=0.1, AL(2018)",".")
ax[:plot](db,convert(Array{Float64,1},sn_big),label="r=0.1, BigFloat",linestyle="--",lw=3)
#ax[:set_xlabel](L"$\log_{10}(1+r-b)$")
ax[:set_title](L"Area of overlap")
ax[:legend](loc="lower right",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-22,1e-2])
x1 = -12.0+cos.(t1); y1=-5.0+aspect*sin.(t1)
x2 = -10.9+0.1*cos.(t1); y2=-5.0+0.1*aspect*sin.(t1)
ax[:plot](x1,y1)
ax[:plot](x2,y2)

ax = axes[4]
ax[:plot](db,abs.(area-convert(Array{Float64,1},sn_big)),label="r=0.1, MA(2002)",lw=3)
ax[:plot](db,abs.(sn-convert(Array{Float64,1},sn_big)),label="r=0.1, AL(2018)",".")
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:set_xlabel](L"$1+r-b$")
#ax[:set_ylabel]("Log Abs(Error in area of overlap)")
ax[:legend](loc="center left",fontsize=8)
ax[:axis]([1e-15,1e-2,1e-25,1e-5])
savefig("area_of_overlap_r01.pdf", bbox_inches="tight")
