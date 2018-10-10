include("occultquad.jl")
include("../../../src/transit_poly_struct.jl")
include("../../../test/loglinspace.jl")

if VERSION >= v"0.7"
  using Statistics
end

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
if b <= abs(R-r)
  area = pi*minimum([r,R])^2
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(abs(sqarea_triangle(R,b,r)))
  # Angle of section for occultor:
  kap0  = atan(kite_area2,(r-R)*(r+R)+b^2)
  # Angle of section for source:
  kap1 = atan(kite_area2,(R-r)*(r+R)+b^2)
  # Flux of visible uniform disk:
  area = R^2*kap1 + r^2*kap0 - 0.5*kite_area2
end
return area
end

# Compute non-linear limb-darkened light curve as a sum over uniform
# sources:

function occult_nonlinear(r::T,b::T,c1::T,c2::T,c3::T,c4::T) where {T <: Real}
#  Computes non-linear limb-darkening with numerical integration
#  using the same approach as Batman:

  function Iofmu(r::T) where {T <: Real}
  mu = sqrt(one(T)-r^2)
  muh = sqrt(mu)
  return one(T)-c1-c2-c3-c4+c1*muh+c2*mu+c3*muh*mu+c4*mu^2
  end

  x0 = maximum([b-r,zero(T)]); xn = minimum([b+r,one(T)])
  theta1 = asin(x0); theta2 = asin(xn);
  nint = 16; delta0 = one(T); delta1 = zero(T)
  tol = 1e-8; iter = 0; itmax = 20
  while abs(delta0-delta1) > tol*delta0 && iter <= itmax
    if VERSION >= v"0.7"
      xgrid = sin.(range(theta1,stop=theta2,length=nint))
    else
      xgrid = sin.(linspace(theta1,theta2,nint))
    end
    delta0 = delta1
    delta1 = zero(T)
    for i=2:nint
      sb = Iofmu(0.5*(xgrid[i]+xgrid[i-1]))
#      println(xgrid[i]," ",xgrid[i-1])
      da = area_overlap(xgrid[i],r,b)-area_overlap(xgrid[i-1],r,b)
      delta1 += sb*da
    end
    nint *=2; iter +=1
#    println("nint: ",nint," delta: ",delta1)
  end
return delta1
end

function transit_poly_int(r,b,c,ns)
if b < (1.0+r)
  s_1 = maximum([0.0,b-r]); s_2=minimum([1.0,b+r])
  ds=(s_2-s_1)/ns
  s = linearspace(s_1+.5*ds,s_2-.5*ds,ns)
  fobs = zero(r)
  for j=1:ns
    if s[j] < r-b
      dphi = 2pi
    else
      dphi = 2*acos((s[j]^2+b^2-r^2)/(2*b*s[j]))
    end
    z = sqrt(sqrt(1-s[j]^2))
    imu = 1.0-sum(c)
#    one(T)-c1-c2-c3-c4+c1*muh+c2*mu+c3*muh*mu+c4*mu^2
    for i=1:4
      imu += c[i]*z^i
    end
    fobs += s[j]*ds*dphi*imu
  end
  norm = 1.-sum(c)
 #omega = (1.0-c1-c2-c3-c4)+0.8*c1+2.0*c2/3.0+4.0*c3/7.0+c4/2.0
  for i=1:4
    norm += 4c[i]/(4.0+i)
  end
  fobs /= norm
#  fobs = 1-fobs
  return fobs
else
  return 1.0
end
end



nx = 2048
if VERSION >= v"0.7"
 x = range(-1.2,stop=1.2,length=nx)
else
 x = linspace(-1.2,1.2,nx)
end
y = 0.0
b = sqrt.(x.^2 .+ y.^2)
c1 = 0.2; c2 = 0.2; c3 = 0.2; c4 = 0.2
r=0.1; flux1 = zeros(nx); flux2 = zeros(nx)
omega = (1.0-c1-c2-c3-c4)+0.8*c1+2.0*c2/3.0+4.0*c3/7.0+c4/2.0
ns = 5000
for i=1:nx
  if b[i] < 1.0+r
    flux1[i]=occult_nonlinear(r,b[i],c1,c2,c3,c4)/omega
    flux2[i]=transit_poly_int(r,b[i],[c1,c2,c3,c4],ns)
#    println(i," r: ",r," b: ",b[i]," flux: ",flux[i])
  end
end

using PyPlot
clf()
plot(x, 1.0 .- flux1/pi)
flux1 = 1.0 .- flux1/pi
plot(x, 1.0 .- flux2/pi)
flux2 = 1.0 .- flux2/pi
# The following two lines are used to check occult_nonlinear
# when c1 = c3 = 0, which yields the quadratic case:
#u2 = -c4; u1 = c2-2*u2
#fquad = occultquad.(b,u1,u2,r)
#plot(x,fquad,linestyle="--")
savefig("occult_nonlinear.pdf", bbox_inches="tight")

# Now, optimize a fit with polynomial limb-darkening:

function optimize_fit_linear!(r::T,b::Array{T,1},fobs::Array{T,1},fmod::Array{T,1},u_n::Array{T,1}) where {T <: Real}
  nb = length(b)   # Number of points in "light curve"
  nu = length(u_n) # Number of limb-darkening coefficients
  # Optimizes a polynomial limb-darkening model to "data":
  # Compute a limb-darkening model with the impact parameter grid b
  # Initialize:
  trans = transit_init(r,b[1],u_n,true)
  # Compute model & gradient at each point:
  smat = zeros(T,nu+1,nu+1); fdots = zeros(T,nu+1)
  for i=1:nb
    trans.b = b[i]
    fcur = transit_poly!(trans)
    for j=1:nu+1
      fdots[j] += fobs[i]*trans.sn[j]
      for k=1:nu+1
        smat[j,k] += trans.sn[j]*trans.sn[k]
      end
    end
  end
  # Now, invert to find coefficients of:  f_obs = \sum_j \alpha_j s_j:
#  println("f.s: ",fdots," s^T s: ",smat)
  alpha = \(smat,fdots)
#  println("alpha: ",alpha)
  fill!(fmod,zero(T))
  for i=1:nb
    trans.b = b[i]
    fcur = transit_poly!(trans)
    for j=1:nu+1
      fmod[i] += alpha[j]*trans.sn[j]
    end
#    println("b: ",b[i]," fobs: ",fobs[i]," fmod: ",fmod[i])
  end
return fmod
end

function optimize_fit_linear!(r::T,b::Array{T,1},fobs::Array{T,1},fmod::Array{T,1},u_n::Array{T,1}) where {T <: Real}
  nb = length(b)   # Number of points in "light curve"
  nu = length(u_n) # Number of limb-darkening coefficients
  # Optimizes a polynomial limb-darkening model to "data":
  # Compute a limb-darkening model with the impact parameter grid b
  # Initialize:
  trans = transit_init(r,b[1],u_n,true)
  # Compute model & gradient at each point:
  smat = zeros(T,nu+1,nu+1); fdots = zeros(T,nu+1)
  for i=1:nb
    trans.b = b[i]
    fcur = transit_poly!(trans)
    for j=1:nu+1
      fdots[j] += fobs[i]*trans.sn[j]
      for k=1:nu+1
        smat[j,k] += trans.sn[j]*trans.sn[k]
      end
    end
  end
  # Now, invert to find coefficients of:  f_obs = \sum_j \alpha_j s_j:
#  println("f.s: ",fdots," s^T s: ",smat)
  alpha = \(smat,fdots)
#  println("alpha: ",alpha)
  fill!(fmod,zero(T))
  for i=1:nb
    trans.b = b[i]
    fcur = transit_poly!(trans)
    for j=1:nu+1
      fmod[i] += alpha[j]*trans.sn[j]
    end
#    println("b: ",b[i]," fobs: ",fobs[i]," fmod: ",fmod[i])
  end
return fmod
end

clf()
#plot(b,flux)
fig,axes = subplots(2,1, figsize=(8,8))
ax = axes[2]
dev = zeros(6)
flux = flux2
fmod = copy(flux)
#u_n = [.1]
#fmod = optimize_fit!(r,b,flux,fmod,u_n)
#plot(b,flux-fmod,linestyle="--",label="linear")
#println("linear: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="quadratic")
println("quad: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="cubic")
println("cubic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1,0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="quartic")
println("quartic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1,0.1,0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="quintic")
println("quintic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1,0.1,0.1,0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="sextic")
println("sextic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
u_n = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]
fmod = optimize_fit!(r,b,flux,fmod,u_n)
ax[:plot](b,(flux-fmod)*1e6,linestyle="--",label="septic")
println("septic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
#u_n = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
#fmod = optimize_fit!(r,b,flux,fmod,u_n)
#plot(b,flux-fmod,linestyle="--",label="octic")
#println("octic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
#u_n = [.1,.1,.1,.1,.1,.1,.1,.1,.1]
#fmod = optimize_fit!(r,b,flux,fmod,u_n)
#plot(b,flux-fmod,linestyle="--",label="nonic")
#println("nonic: ",maximum(abs,flux-fmod)," sig: ",std(flux-fmod))
ax[:legend](loc = "upper left",fontsize=8)
ax[:set_ylabel]("Deviation of fit [ppm]")
ax[:set_xlabel](L"$b$", fontsize=18)
ax[:axis]([0,1.2,-10,10])
ax = axes[1]

ax[:plot](b,flux,linewidth=2,label="non-linear")
ax[:plot](b,fmod,linestyle="--",label="septic",linewidth=2)
ax[:set_title](L"$c_1 = c_2=c_3=c_4=0.2$")
ax[:legend](loc="upper left", fontsize=8)
ax[:axis]([0,1.2,0.984,1.001])
ax[:set_ylabel]("Relative flux")

show()

savefig("occult_nonlinear_poly.pdf",bbox_inches="tight")
