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
sn_big = convert(Float64,kap1_big + r_big^2*kap0_big - .5*kite_area2_big)

return area,sn,sn_big
end


using PyPlot
fig,axes = subplots(2,2)

r = 0.1
nb = 1000
db = logspace(-15,-2,nb)
b = 1-r+db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[1]
ax[:plot](log10.(db),log10.(abs.(pi*r^2-area)),label="r=0.01, MA(2002)",lw=3)
ax[:plot](log10.(db),log10.(abs.(pi*r^2-sn)),label="r=0.01, RA(2018)",linestyle=":",lw=3)
ax[:plot](log10.(db),log10.(abs.(pi*r^2-sn_big)),label="r=0.01, BigFloat",linestyle="--")
#ax[:set_xlabel](L"$\log_{10}(b-(1-r))$")
ax[:set_title ](L"$\log_{10}(\pi r^2$-Area of overlap)")
ax[:legend](loc="lower right",fontsize=8)

ax = axes[2]
ax[:plot](log10.(db),log10.(abs.(area-sn_big)),label="r=0.01, MA(2002)",lw=3)
ax[:plot](log10.(db),log10.(abs.(sn-sn_big)),label="r=0.01, RA(2018)",".")
ax[:set_xlabel](L"$\log_{10}(b-(1-r))$")
ax[:set_ylabel]("Log Abs(Error in area of overlap)")
ax[:legend](loc="left",fontsize=8)

b = 1+r-db
area=zeros(nb); sn=zeros(nb); sn_big=zeros(nb)
for i=1:nb
  area[i],sn[i],sn_big[i] = circle_overlap(r,b[i])
end

ax = axes[3]
mask = area .> 0.
ax[:plot](log10.(db[mask]),log10.(area[mask]),label="r=0.01, MA(2002)",lw=3)
ax[:plot](log10.(db),log10.(sn),label="r=0.01, RA(2018)",linestyle=":",lw=3)
ax[:plot](log10.(db),log10.(sn_big),label="r=0.01, BigFloat",linestyle="--")
#ax[:set_xlabel](L"$\log_{10}(1+r-b)$")
ax[:set_title]("Log Area of overlap")
ax[:legend](loc="lower right",fontsize=8)

ax = axes[4]
ax[:plot](log10.(db),log10.(abs.(area-sn_big)),label="r=0.01, MA(2002)",lw=3)
ax[:plot](log10.(db),log10.(abs.(sn-sn_big)),label="r=0.01, RA(2018)",".")
ax[:set_xlabel](L"$\log_{10}(1+r-b)$")
#ax[:set_ylabel]("Log Abs(Error in area of overlap)")
ax[:legend](loc="left",fontsize=8)

savefig("area_of_overlap_r0.1.pdf", bbox_inches="tight")
