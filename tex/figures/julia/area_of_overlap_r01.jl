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
fig,axes = subplots(2,2,figsize=(12, 7))
fig[:subplots_adjust](wspace=0.25)

r = 0.1
nb = 1000
db = logspace(-15,-2,nb)
b = 1-r+db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(BigFloat,nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[1]
ax[:plot](db,abs.(pi*r^2-area),label=L"$r=0.1$, MA(2002)",lw=2)
ax[:plot](db,abs.(pi*r^2-sn),label=L"$r=0.1$, AL(2018)",".")
ax[:plot](db,abs.(convert(Array{Float64,1},pi*big(r)^2-sn_big)),label=L"$r=0.1$, BigFloat",linestyle="--",lw=2)
ax[:set_ylabel](L"$\pi r^2 - A_\mathrm{overlap}$", fontsize=14)
ax[:legend](loc="lower right",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-25,1e-2])
ax[:plot]([1e-16,1e-2],[1,1]*(pi*0.1^2/2^53),c="grey",linestyle="-.")

ax = fig[:add_axes]((0.25, 0.925, 0.1, 0.1))
t1 = linspace(0,2pi,100)
x1 = cos.(t1); y1 = sin.(t1)
x2 = 0.8 + 0.2 * cos.(t1); y2 = 0.2 * sin.(t1)
ax[:plot](x1,y1,"k-",lw=1.5)
ax[:fill_between](x2,0,y2,facecolor="k")
ax[:set_aspect](1)
ax[:axis]("off")

ax = fig[:add_axes]((0.685, 0.925, 0.1, 0.1))
t1 = linspace(0,2pi,100)
x1 = cos.(t1); y1 = sin.(t1)
x2 = 1.2 + 0.2 * cos.(t1); y2 = 0.2 * sin.(t1)
ax[:plot](x1,y1,"k-",lw=1.5)
ax[:fill_between](x2,0,y2,facecolor="k")
ax[:set_aspect](1)
ax[:axis]("off")

ax = axes[2]
ax[:plot](db,abs.(area-convert(Array{Float64,1},sn_big)),label=L"$r=0.1$, MA(2002)",lw=2)
ax[:plot](db,abs.(sn-convert(Array{Float64,1},sn_big)),label=L"$r=0.1$, AL(2018)",".")
ax[:set_xlabel](L"$b-(1-r)$", fontsize=14)
ax[:set_ylabel]("error", fontsize=14)
ax[:legend](loc="upper right",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-25,1e-2])
ax[:plot]([1e-16,1e-2],[1,1]*(pi*0.1^2/2^53),c="grey",linestyle="-.")

b = 1+r-db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(BigFloat,nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[3]
mask = area .> 0.0
ax[:plot](db[mask],area[mask],label=L"$r=0.1$, MA(2002)",lw=2)
ax[:plot](db,sn,label=L"$r=0.1$, AL(2018)",".")
ax[:plot](db,convert(Array{Float64,1},sn_big),label=L"$r=0.1$, BigFloat",linestyle="--",lw=2)
ax[:set_ylabel](L"$A_\mathrm{overlap}$", fontsize=14)
ax[:legend](loc="lower right",fontsize=8)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:axis]([1e-15,1e-2,1e-25,1e-2])

ax = axes[4]
ax[:plot](db,abs.(area-convert(Array{Float64,1},sn_big)),label=L"$r=0.1$, MA(2002)",lw=2)
ax[:plot](db,abs.(sn-convert(Array{Float64,1},sn_big)),label=L"$r=0.1$, AL(2018)",".")
ax[:set_ylabel]("error", fontsize=14)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:set_xlabel](L"$1+r-b$", fontsize=14)
ax[:legend](loc="upper right",fontsize=8)
ax[:axis]([1e-15,1e-2,1e-25,1e-2])
savefig("area_of_overlap_r01.pdf", bbox_inches="tight")
