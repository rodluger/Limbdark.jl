include("occultquad.jl")
include("../../../src/transit_poly_struct.jl")

function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}
# How to compute (sixteen times the) area squared of triangle with
# high precision (Goldberg 1991).
# First, do a quick sort of three numbers:
if c > b
  tmp = c; c = b; b = tmp
end
if c > a
  tmp = a; a = c; c = tmp
end
if b > a
  tmp = b; b = a; a = tmp
end
area = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
return area
end

function area_overlap(R::T,r::T,b::T) where {T <: Real}
# Compute the light curve of a uniform source:
if b <= R-r
  area = pi*r^2
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(abs(sqarea_triangle(R,b,r)))
  # Angle of section for occultor:
  kap0  = atan2(kite_area2,(r-R)*(r+R)+b^2)
  # Angle of section for source:
  kap1 = atan2(kite_area2,(R-r)*(r+R)+b^2)
  # Flux of visible uniform disk:
  area = R^2*kap1 + r^2*kap0 - .5*kite_area2
end
return area 
end

# Compute non-linear limb-darkened light curve as a sum over uniform
# sources:

function occult_nonlinear(r::T,b::T,c1::T,c2::T,c3::T,c4::T) where {T <: Real}
#  Computes non-linear limb-darkening with numerical integration:

  function Iofmu(r::T) where {T <: Real}
  mu = sqrt(one(T)-r^2)
  muh = sqrt(mu)
  return one(T)-c1-c2-c3-c4+c1*muh+c2*mu+c3*muh*mu+c4*mu^2
  end
  
  x0 = maximum([b-r,zero(T)]); xn = minimum([b+r,one(T)])
  theta1 = asin(x0); theta2 = asin(xn);
  nint = 16; delta0 = one(T); delta1 = zero(T)
  tol = 1e-6; iter = 0; itmax = 20
  while abs(delta0-delta1) > tol*delta0 && iter <= itmax
    xgrid = sin.(linspace(theta1,theta2,nint))
    delta0 = delta1
    delta1 = zero(T)
    for i=2:nint
      sb = Iofmu(.5*(xgrid[i]+xgrid[i-1]))
#      println(xgrid[i]," ",xgrid[i-1])
      da = area_overlap(xgrid[i],r,b)-area_overlap(xgrid[i-1],r,b)
      delta1 += sb*da
    end
    nint *=2; iter +=1
#    println("nint: ",nint," delta: ",delta1)
  end
return delta1
end

nx = 1024
x = linspace(-1.2,1.2,nx)
y = 0.5
b = sqrt.(x.^2+y.^2)
c1 = 0.2; c2 = 0.2; c3 = 0.2; c4 = 0.2
r=0.1; flux = zeros(nx)
omega = (1.-c1-c2-c3-c4)+.8*c1+2.*c2/3.+4.*c3/7.+c4/2.
for i=1:nx
  if b[i] < 1.0+r
    flux[i]=occult_nonlinear(r,b[i],c1,c2,c3,c4)/omega
#    println(i," r: ",r," b: ",b[i]," flux: ",flux[i])
  end
end

using PyPlot
clf()
plot(x,1.-flux/pi)
# The following two lines are used to check occult_nonlinear
# when c1 = c3 = 0, which yields the quadratic case:
#u2 = -c4; u1 = c2-2*u2
#fquad = occultquad.(b,u1,u2,r)
#plot(x,fquad,linestyle="--")
savefig("occult_nonlinear.pdf", bbox_inches="tight")

# Now, optimize a fit with polynomial limb-darkening:

function optimize_fit!(r::T,b::Array{T,1},fobs::Array{T,1},fmod::Array{T,1},u_n::Array{T,1}) where {T <: Real}
  nb = length(b)   # Number of points in "light curve"
  nu = length(u_n) # Number of limb-darkening coefficients
  # Optimizes a polynomial limb-darkening model to "data":
  # Compute a limb-darkening model with the impact parameter grid b
  # Initialize:
  trans = transit_init(r,b[1],u_n,true)
  # Compute model & gradient at each point:
  smat = zeros(T,nu+1,nu+1); fdots = zeros(T,nu+1)
  for i=1:b
    trans.b = b[i]
    flux = transit_poly!(trans)
    for j=1:nu+1
      fdots[j] += fobs[i]*trans.sn[j]
      for k=1:nu+1
        smat[j,k] += trans.sn[j]*trans.sn[k]
      end
    end
  end
  # Now, invert to find coefficients of:  f_obs = \sum_j \alpha_j s_j:
  alpha = \(smat,fdots)
return
end
