# Define constants for speed:
include("define_constants.jl")
# Include linear algebara:

if VERSION >= v"0.7"
  using LinearAlgebra
end

# Include definition of Transit structure type:
include("transit_structure.jl")
# Include code which computes linear limb-darkening term:
include("s2.jl")
# Include code which computes M_m:
include("Mm_compute.jl")

# Computes the coefficient for the uniform disk case, S_0:
function compute_uniform!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; r2=r*r; b2=b*b
# Uniform disk case:
# Compute sn[1] and its derivatives:
if b <= 1-r  # k^2 > 1
  t.sn[1] = pi*(1-r2)
  if t.grad
    t.dsndr[1] = -2*pi*r
    t.dsndb[1] = 0.
  end
  t.kap0 = convert(T,pi); t.kck = zero(T)
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
#  t.kite_area2 = sqrt(sqarea_triangle(one(T),b,r)) 
  t.kite_area2 = sqrt(t.sqarea)
  # Angle of section for occultor:
  t.kap0  = atan(t.kite_area2,(r-1)*(r+1)+b2)
  # Angle of section for source:
  t.pimkap1 = atan(t.kite_area2,(r-1)*(r+1)-b2)
  # Flux of visible uniform disk:
  t.sn[1] = t.pimkap1 - r2*t.kap0 + .5*t.kite_area2
  if t.grad
    t.dsndr[1]= -2*r*t.kap0
    t.dsndb[1]= t.kite_area2/b
  end
  t.kck = t.kite_area2*t.fourbrinv
end
return
end

function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}
# Function which computes sixteen times the square of the area
# of a triangle with sides a, b and c using Kahan method.
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

# Computes a limb-darkened transit light curve with the dependence:
# I(\mu) = 1-\sum_{n=1}^N u_n (1-\mu)^n
# where \mu = \cos{\theta} = z is the cosine of the angle
# angle from the sub-stellar point (or equivalently the
# height on the star relative to the sky plane if the radius 
# of the star is unity.

function transit_poly_c(t::Transit_Struct{T}) where {T <: Real}
# Number of limb-darkening components to include (beyond 0 and 1):
# We are parameterizing these with the function:
# g_n = c_n [(n+2) z^n - n z^{n-2}] for n >= 2
# while g_{0} = c_0 z^0 (uniform source) and g_{1} = c_1 z^1 (linear limb-darkening)
# which gives a Green's integral of:
# P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
# for which we have a solution in terms of M_m.
# The variable "t" is a structure which contains transit parameters
# and intermediate quantities computed from these:
r=t.r; b=t.b; n = t.n; r2=r*r; b2 = b*b
# Set up a vector for storing results of P(G_n)-Q(G_n); note that
# this is a different vector than the Starry case:

# Check for different cases:
if b >= 1+r
  # unobscured - return one:
  return one(T)
end
if r >= 1+b
  # full obscuration - return zero:
  return zero(T)
end
if b == 0.0
  # Annular eclipse - integrate around the full boundary of both bodies:
  onemr2 = 1-r2
  flux = zero(T); t.sqrt1mr2 = sqrt(onemr2)
  flux = (t.c_n[1]*onemr2+t.twothird*t.c_n[2]*t.sqrt1mr2^3)
  fac= 2r2*onemr2
  @inbounds for i=2:t.n
    flux += -t.c_n[i+1]*fac
    fac *= t.sqrt1mr2
  end
  return flux*pi*t.den
else
# Next, compute k^2 = m:
  t.onembmr2=(r+1-b)*(1-r+b); t.fourbr = 4b*r; t.fourbrinv = inv(t.fourbr)
  t.sqbr = sqrt(b*r); t.sqbrinv = inv(t.sqbr)
  t.onembmr2inv=inv(t.onembmr2); t.sqonembmr2 = sqrt(t.onembmr2)
  t.onembpr2 = (1-r-b)*(1+b+r)
  t.sqarea = sqarea_triangle(one(T),r,b)
#  t.k2 = t.onembmr2*t.fourbrinv
  t.k2 = t.onembpr2*t.fourbrinv+1
  t.onemr2mb2 = (1.0-r)*(1.0+r)-b2
  if t.k2 > 0
    t.k = sqrt(t.k2)
  else
    println("negative k2: ",t.k2," r: ",r," b: ",b)
    t.k2 = 0.0; t.k = 0.0
  end
  if t.k2 > 1
    if t.k2 > 2.0
      t.kc2 = 1.0-inv(t.k2)
      t.kc = sqrt(t.kc2)
    else
      t.kc2 = t.onembpr2/((1-b+r)*(1-r+b))
      t.kc = sqrt(t.kc2)
    end
  else
    if t.k2 > 0.5
      t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
      t.kc = sqrt(t.kc2)
    else
      t.kc2 = 1.0-t.k2
      t.kc = sqrt(t.kc2)
    end
  end
end

# Compute uniform case, sn[1]:
compute_uniform!(t)
# If uniform model, then return:
if t.n == 0
  return t.c_n[1]*t.sn[1]*t.den
end

# Compute linear case, sn[2]:
#s2!(t)
t.sn[2],t.Eofk,t.Em1mKdm = s2_ell(r,b)
#sn2,Eofk,Em1mKdm = s2_ell(r,b)
# Now, compare results:
#if abs(t.sn[2]-sn2) > 1e-8*abs(sn2)
#  println("r: ",r," b: ",b," sn[2]: ",t.sn[2]," sn2: ",sn2)
#end
# Special case of linear limb-darkening:
if t.n == 1
  flux = t.c_n[1]*t.sn[1]+t.c_n[2]*t.sn[2]
  flux *= t.den  # for c_2 and above, the flux is zero.
  return flux
end

# Special case of quadratic limb-darkening:
if t.n == 2
# Transformed expressions from Mandel & Agol:
  eta2 = 0.5*r2*(r2+2*b2)
  if t.k2 > 1
    four_pi_eta = 4pi*eta2
  else
    four_pi_eta = 2*(pi-t.pimkap1+2*eta2*t.kap0-0.25*t.kite_area2*(1.0+5r2+b2))
  end
  t.sn[3] = 2*(t.sn[1]-pi)+four_pi_eta
  flux = t.c_n[1]*t.sn[1]+t.c_n[2]*t.sn[2]+t.c_n[3]*t.sn[3]
  flux *= t.den  # for c_2 and above, the flux is zero.
  return flux
end

# Compute the M_m functions:
if t.k2 > 0
#  if (t.k2 < 0.5 || t.k2 > 250.0) && t.n > 3
  if (t.k2 < 0.5) && t.n > 3
# This computes M_m for the largest four m, and then works down to smaller values:
    Mm_lower!(t)
  else
# This computes Mm for m=0 to 3, and then works upward to larger m:
    Mm_raise!(t)
  end
end

# Add up first two terms in flux numerator:
flux = t.c_n[1]*t.sn[1]+t.c_n[2]*t.sn[2]
# Next, loop over the Green's function components:
@inbounds for n=2:t.n
#  pofgn_M = (1+(r-b)*(r+b))*t.Mm[n+1]-t.Mm[n+3]
#  pofgn_M = t.onemr2mb2*t.Mm[n+1]-t.Mm[n+3]
#  pofgn_M = 2*r^2*t.Mm[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mm[n+1]+(1-(b-r)^2)*((b+r)^2-1)*t.Mm[n-1])
#  pofgn_M = 2*r^2*t.Mm[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mm[n+1]+sqarea_triangle(one(T),r,b)*t.Mm[n-1])
#  pofgn_M = 2*r^2*t.Mm[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mm[n+1]+t.sqarea*t.Mm[n-1])
#  pofgn_M = 2*r2*t.Mm[n+1]-n/(n+2)*(t.onemr2mb2*t.Mm[n+1]+t.sqarea*t.Mm[n-1])
  pofgn_M = 2*r2*t.Mm[n+1]-n*t.minv[n+2]*(t.onemr2mb2*t.Mm[n+1]+t.sqarea*t.Mm[n-1])
# Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
# boundary for n > 0.
# Compute sn[n]:
  t.sn[n+1] = -pofgn_M
  if t.b <= 1e-6 && r < 1
    # Use analytic formula near b=0:
    t.sqrt1mr2 = sqrt(1-r2)
    t.sn[n+1] = -0.5*pi*r2*t.sqrt1mr2^(n-4)*(4*(1-r2)^2+n*b2*((2+n)*r2-4))
  end
  flux += t.c_n[n+1]*t.sn[n+1]
end
# That's it!
# flux = sum(t.c_n.*t.sn)/(pi*(t.c_n[1]+t.twothird*t.c_n[2]))  # for c_2 and above, the flux is zero.
# Divide by denominator:
flux *= t.den  # for c_2 and above, the flux is zero.
return flux
end
# That's it!

function transit_poly!(r::T,b::T,u_n::Array{T,1},dfdrbu::Array{T,1}) where {T <: Real}
t = transit_init(r,b,u_n,true)
# Pass c_n (without last two dummy values):
flux = transit_poly_c!(t)
# Now, transform derivaties from c to u:
fill!(dfdrbu,zero(T))
dfdrbu[1] = t.dfdrbc[1]  # r derivative
dfdrbu[2] = t.dfdrbc[2]  # b derivative
# u_n derivatives:
@inbounds for i=1:t.n, j=0:t.n
  dfdrbu[i+2] += t.dfdrbc[j+3]*t.dcdu[j+1,i]
end
return flux
end

function transit_poly(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
t = transit_init(r,b,u_n,false)
# Pass c_n (without last two dummy values):
return transit_poly_c(t) 
end

function transit_poly!(t::Transit_Struct{T}) where {T <: Real}
# This function assumes that c_n has already been computed from u_n
# (this can be used to save compute time when limb-darkening is fixed for
# a range of radii/impact parameters).
# Pass transit structure, and compute flux:
if t.grad
  flux = transit_poly_c!(t)
  # Now, transform derivaties from c to u:
  fill!(t.dfdrbu,zero(T))
  t.dfdrbu[1] = t.dfdrbc[1]  # r derivative
  t.dfdrbu[2] = t.dfdrbc[2]  # b derivative
#  t.dfdrbu[3:t.n+2]=BLAS.gemv!('T',1.0,t.dcdu,t.dfdrbc[3:t.n+3],0.0,t.dfdrbu[3:t.n+2])
  BLAS.gemv!('T',1.0,t.dcdu,t.dfdrbc[3:t.n+3],0.0,t.dfdrbu[3:t.n+2])
#  t.dfdrbu[3:t.n+2]=BLAS.gemv('T',1.0,t.dcdu,t.dfdrbc[3:t.n+3])
#  @inbounds for i=1:t.n, j=0:t.n
#    t.dfdrbu[i+2] += t.dfdrbc[j+3]*t.dcdu[j+1,i]
#  end
  return flux
else
  return transit_poly_c(t)
end
return
end

function transit_poly_c!(t::Transit_Struct{T}) where {T <: Real}
r = t.r; b=t.b; n = t.n; r2 =r*r; b2=b*b
@assert((length(t.c_n)+2) == length(t.dfdrbc))
@assert(r > 0)
# Number of limb-darkening components to include (beyond 0 and 1):
# We are parameterizing these with the function:
# g_n = c_n [(n+2) z^n - n z^{n-2}] for n >= 2
# while g_{0} = c_0 z^0 (uniform source) and g_{1} = c_1 z^1 (linear limb-darkening)
# which gives a Green's integral of:
# P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
# for which we have a solution in terms of M_m.
# Compute the derivative of the flux with respect to the different coefficients.

# Set up a vector for storing results of P(G_n)-Q(G_n); note that
# this is a different vector than the Starry case:
# Check for different cases:
if b >= 1+r || r ==  0.0
  # unobscured - return one, and zero derivatives:
  fill!(t.dfdrbc,zero(T))
  return one(T)
end
if r >= 1+b
  # full obscuration - return zero, and zero derivatives:
  fill!(t.dfdrbc,zero(T))
  return zero(T)
end
if b == 0.0
  # Annular eclipse - integrate around the full boundary of both bodies:
  flux = zero(T); onemr2 = 1-r2; t.sqrt1mr2 = sqrt(onemr2)
  fill!(t.dfdrbc,zero(T))
  flux = (t.c_n[1]*onemr2+t.twothird*t.c_n[2]*t.sqrt1mr2^3)*pi*t.den
  fac  = 2r2*onemr2*pi*t.den
  facd = -2r*pi*t.den
  t.dfdrbc[1] = t.c_n[1]*facd + t.c_n[2]*facd*t.sqrt1mr2
  @inbounds for i=2:t.n
    flux -= t.c_n[i+1]*fac
    t.dfdrbc[1] += t.c_n[i+1]*facd*(2*onemr2-i*r2)
    t.dfdrbc[i+3] -= fac
    fac *= t.sqrt1mr2
    facd *= t.sqrt1mr2
  end
  #  dfdrbc[2]=0 since the derivative with respect to b is zero.
  t.dfdrbc[3] = (onemr2-flux)*pi*t.den
  t.dfdrbc[4] = t.twothird*(t.sqrt1mr2^3-flux)*pi*t.den
  # Also need to compute derivatives [ ]
  return flux
else
# Next, compute k^2 = m:
  t.onembmr2=(r-b+1)*(1-r+b); t.fourbr = 4b*r; t.fourbrinv = inv(t.fourbr)
  t.sqbr = sqrt(b*r); t.sqbrinv = inv(t.sqbr)
  t.onembmr2inv=inv(t.onembmr2); t.sqonembmr2 = sqrt(t.onembmr2)
  t.onembpr2 = (1-r-b)*(1+b+r)
  t.sqarea = sqarea_triangle(one(T),r,b)
  t.k2 = t.onembmr2*t.fourbrinv; 
  t.onemr2mb2 = (1.0-r)*(1.0+r)-b2
  if t.k2 > 0
    t.k = sqrt(t.k2)
  else
    println("negative k2: ",t.k2," r: ",r," b: ",b)
    t.k2 = 0.0; t.k = 0.0
  end
  dkdr = (b2-r2-1)/(8*t.k*b*r2)
  dkdb = (r2-b2-1)/(8*t.k*b2*r)
  if t.k2 > 1
    if t.k2 > 2.0
      t.kc2 = 1.0-inv(t.k2)
      t.kc = sqrt(t.kc2)
    else
      t.kc2 = t.onembpr2/((1+r-b)*(1-r+b))
      t.kc = sqrt(t.kc2)
    end
  else
    if t.k2 > 0.5
      t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
      t.kc = sqrt(t.kc2)
    else
      t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
      t.kc = sqrt(t.kc2)
    end
  end
end

# Compute uniform case:
compute_uniform!(t)

# Compute sn[2] and its derivatives:
#s2!(t)
if t.n >= 1
  t.sn[2],t.Eofk,t.Em1mKdm = s2!(r,b,t.s2_grad)
  t.dsndr[2] = t.s2_grad[1]
  t.dsndb[2] = t.s2_grad[2]
end


# Special case of quadratic limb-darkening:
if t.n >= 2
  if t.n == 2
# Transformed expressions from Mandel & Agol:
    r2pb2 = r2+b2
    eta2 = 0.5*r2*(r2pb2+b2)
    deta2dr =  2*r*r2pb2
    deta2db = 2*b*r2
    if t.k2 > 1
      four_pi_eta = 4pi*(eta2-0.5)
      detadr = 4pi*deta2dr
      detadb = 4pi*deta2db
    else
      four_pi_eta = 2*(-t.pimkap1+2*eta2*t.kap0-0.25*t.kite_area2*(1.0+5r2+b2))
      detadr = 8r*(r2pb2*t.kap0-t.kite_area2)
      detadb = 2/b*(4*b2*r2*t.kap0-(1+r2pb2)*t.kite_area2)
    end
    t.sn[3] = 2*t.sn[1]+four_pi_eta
    t.dsndr[3] = 2*t.dsndr[1]+detadr
    t.dsndb[3] = 2*t.dsndb[1]+detadb
  else
    # Compute the M_m functions:
    if t.k2 > 0
#      if (t.k2 < 0.5 || t.k2 > 2.0) && t.v_max > 3
      if (t.k2 < 0.5) && t.n > 3
    # This computes Mm for largest four m, and then works down to smaller values:
        Mm_lower!(t)
      else
    # This computes Mm and then works upward to larger m:
        Mm_raise!(t)
      end
    end
  
    # Next, loop over the Green's function components:
    binv = inv(b)
    @inbounds for n=2:t.n
#      pofgn_M = (1+(r-b)*(r+b))*t.Mm[n+1]-t.Mm[n+3]
#      pofgn_M = t.onemr2mb2*t.Mm[n+1]-t.Mm[n+3]
#      pofgn_M = 2*r^2*t.Mm[n+1]-n/(n+2)*((1-r2-b2)*t.Mm[n+1]+t.kite_area2^2*t.Mm[n-1])
#      pofgn_M = 2*r^2*t.Mm[n+1]-n/(n+2)*(t.onemr2mb2*t.Mm[n+1]+sqarea_triangle(one(T),r,b)*t.Mm[n-1])
#      pofgn_M = 2*r2*t.Mm[n+1]-n/(n+2)*(t.onemr2mb2*t.Mm[n+1]+t.sqarea*t.Mm[n-1])
      pofgn_M = 2*r2*t.Mm[n+1]-n*t.minv[n+2]*(t.onemr2mb2*t.Mm[n+1]+t.sqarea*t.Mm[n-1])
    # Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
    # boundary for n > 0.
    # Compute sn[n]:
      t.sn[n+1] = -pofgn_M
      dpdr_M = 2*r*((n+2)*t.Mm[n+1]-n*t.Mm[n-1])
      t.dsndr[n+1] = -dpdr_M
      dpdb_M = n*binv*((t.Mm[n+1]-t.Mm[n-1])*(r2+b2)+(b2-r2)^2*t.Mm[n-1])
      t.dsndb[n+1] = -dpdb_M
      # When b is very small and r < 1, we'll use a Taylor-series expansion to O(b^3):
      if t.b <= 1e-6 && r < 1
        # Use analytic formula near b=0:
        t.sqrt1mr2 = sqrt(1-r2)
        t.sn[n+1] = -0.5*pi*r2*t.sqrt1mr2^(n-4)*(4*(1-r2)^2+n*b2*((2+n)*r2-4))
        t.dsndr[n+1] = 0.5*pi*r*t.sqrt1mr2^(n-6)*(4*(1-r2)^2*((2+n)*r2-2)+
                         n*b2*(8-8*n*r2+n*(2+n)*r2*r2))
        t.dsndb[n+1] = pi*r2*n*b*t.sqrt1mr2^(n-4)*(4-(2+n)*r2)
      end
      # Handle r=1, n=2 case for small b values:
      if t.b <= 1e-8 && r == 1.0 && n == 2
        t.dsndb[n+1] = -4.0
      end
    end
  end
end
# That's it!
# Compute derivatives with respect to the coefficients:
flux = zero(T)
t.dfdrbc[1]=zero(T)  # Derivative with respect to r
t.dfdrbc[2]=zero(T)  # Derivative with respect to b
@inbounds for i=0:t.n
  # derivatives with respect to the coefficients:
  t.dfdrbc[i+3]= t.sn[i+1]*t.den
  # total flux:
  flux += t.c_n[i+1]*t.dfdrbc[i+3]
  # derivatives with respect to r and b:
  t.dfdrbc[1] += t.c_n[i+1]*t.dsndr[i+1]*t.den
  t.dfdrbc[2] += t.c_n[i+1]*t.dsndb[i+1]*t.den
end
# Include derivatives with respect to first two c_n parameters:
t.dfdrbc[3] -= flux*t.den*pi
t.dfdrbc[4] -= flux*t.den*pi*t.twothird
return flux
end
