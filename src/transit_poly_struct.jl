include("transit_structure.jl")
include("compute_c_n_struct.jl")
include("s2.jl")
include("IJv_derivative_struct.jl")

function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}
# How to compute (twice) area squared of triangle with 
# high precision (Goldberg 1991):
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

function sqarea_triangle(x::Array{T,1}) where {T <: Real}
# How to compute (twice) area squared of triangle with 
# high precision (Goldberg 1991):
a=maximum(x); b=median(x); c=minimum(x)
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
# for which we have a solution in terms of I_v (for even n) and J_v (for odd n).
# The variable "t" is a structure which contains transit parameters
# and intermediate quantities computed from these:
r=t.r; b=t.b; n = t.n
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
  flux = zero(T); sqrt1mr2 = sqrt(1-r^2)
  flux = (t.c_n[1]*(1-r^2)+2/3*t.c_n[2]*sqrt1mr2^3)
  fac= 2r^2*(1-r^2)
  for i=2:t.n
    flux += -t.c_n[i+1]*fac
    fac *= sqrt1mr2
  end
  return flux/(t.c_n[1]+2*t.c_n[2]/3)
else
# Next, compute k^2 = m:
  onembmr2=(r+1-b)*(1-r+b); fourbr = 4b*r
  k2 = onembmr2/fourbr
  if k2 > 1
    if k2 > 2.0
      kc = sqrt(1.-inv(k2))
    else
      kc2 = (1-r-b)*(1+b+r)/(1-b+r)/(1-r+b)
      kc = sqrt(kc2)
    end
  else
    if k2 > 0.5
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
    else
      kc = sqrt(1.-k2)
    end
  end
end

# Compute the highest value of v in J_v or I_v that we need:
# Compute sn[1] and sn[2]:
# Uniform disk case:
if b <= 1-r
  lam = pi*r^2
  t.sn[1] = pi-lam
  kap0 = convert(T,pi); kck = zero(T)
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(sqarea_triangle(one(r),b,r))
  # Angle of section for occultor:
  kap0  = atan2(kite_area2,(r-1)*(r+1)+b^2)
  # Angle of section for source:
  pimkap1 = atan2(kite_area2,(r-1)*(r+1)-b^2)
  # Flux of visible uniform disk:
  t.sn[1] = pimkap1 - r^2*kap0 + .5*kite_area2
  kck = kite_area2/(4*b*r)
end
t.sn[2],Eofk,Em1mKdm = s2_ell(r,b)
# Compute the J_v and I_v functions:
if k2 > 0
  if k2 < 0.5 || k2 > 2.0
# This computes I_v,J_v for the largest v, and then works down to smaller values:
    IJv_lower!(k2,kck,kc,kap0,Eofk,Em1mKdm,t)
  else
# This computes I_0,J_0,J_1, and then works upward to larger v:
    IJv_raise!(k2,kck,kc,kap0,Eofk,Em1mKdm,t)
  end
end

# Next, loop over the Green's function components:
for n=2:t.n
  pofgn = zero(T)
  if iseven(n)
# For even values of n, sum over I_v:
    n0 = convert(Int64,n/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*t.Iv[n0+1]+2b*t.Iv[n0+2])
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
      pofgn += coeff*((r-b)*t.Iv[n0-i+1]+2b*t.Iv[n0-i+2])
    end
    pofgn *= 2r
  else
# Now do the same for odd N_c in sum over J_v:
    n0 = convert(Int64,(n-3)/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*t.Jv[n0+1]+2b*t.Jv[n0+2])
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
      pofgn += coeff*((r-b)*t.Jv[n0-i+1]+2b*t.Jv[n0-i+2])
    end
    pofgn *= 2r*onembmr2^1.5
  end
# Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
# boundary for n > 0.
# Compute sn[n]:
  t.sn[n+1] = -pofgn
end
# That's it!
flux = sum(t.c_n.*t.sn)/(pi*(t.c_n[1]+2*t.c_n[2]/3))  # for c_2 and above, the flux is zero.
return flux
end
# That's it!

function transit_poly!(r::T,b::T,u_n::Array{T,1},dfdrbu::Array{T,1}) where {T <: Real}
t = transit_init(r,b,u_n,true)
#compute_c_n_grad!(t) #This is now included in transit_init function
# Pass c_n (without last two dummy values):
flux = transit_poly_c!(t)
# Now, transform derivaties from c to u:
fill!(dfdrbu,zero(T))
dfdrbu[1] = t.dfdrbc[1]  # r derivative
dfdrbu[2] = t.dfdrbc[2]  # b derivative
# u_n derivatives:
for i=1:t.n, j=0:t.n
  dfdrbu[i+2] += t.dfdrbc[j+3]*t.dcdu[j+1,i]
end
return flux
end

function transit_poly(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
t = transit_init(r,b,u_n,false)
#t.c_n = compute_c_n(t) #This is now included in transit_init function
# Pass c_n (without last two dummy values):
return transit_poly_c!(t) 
end

function transit_poly!(t::Transit_Struct{T}) where {T <: Real}
# This function assumes that c_n has already been computed from u_n
# (this can be used to save compute time when limb-darkening is fixed for
# a range of radii/impact parameters).
# Pass transit structure, and compute flux:
flux = transit_poly_c!(t)
# Now, transform derivaties from c to u:
fill!(t.dfdrbu,zero(T))
t.dfdrbu[1] = t.dfdrbc[1]  # r derivative
t.dfdrbu[2] = t.dfdrbc[2]  # b derivative
# u_n derivatives:
for i=1:t.n, j=0:t.n
  t.dfdrbu[i+2] += t.dfdrbc[j+3]*t.dcdu[j+1,i]
end
return flux
end

function transit_poly_c!(t::Transit_Struct{T}) where {T <: Real}
r = t.r; b=t.b; n = t.n
@assert((length(t.c_n)+2) == length(t.dfdrbc))
@assert(r > 0)
# Number of limb-darkening components to include (beyond 0 and 1):
# We are parameterizing these with the function:
# g_n = c_n [(n+2) z^n - n z^{n-2}] for n >= 2
# while g_{0} = c_0 z^0 (uniform source) and g_{1} = c_1 z^1 (linear limb-darkening)
# which gives a Green's integral of:
# P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
# for which we have a solution in terms of I_v (for even n) and J_v (for odd n).
# Compute the derivative of the flux with respect to the different coefficients.

# Set up a vector for storing results of P(G_n)-Q(G_n); note that
# this is a different vector than the Starry case:
# Check for different cases:
if b >= 1+r || r ==  0.0
  # unobscured - return one:
  return one(T)
end
if r >= 1+b
  # full obscuration - return zero:
  return zero(T)
end
if b == 0.0
  # Annular eclipse - integrate around the full boundary of both bodies:
  flux = zero(T); onemr2 = 1-r^2; sqrt1mr2 = sqrt(onemr2)
  fill!(t.dfdrbc,zero(T))
  den = inv(t.c_n[1]+2*t.c_n[2]/3)
  flux = (t.c_n[1]*onemr2+2/3*t.c_n[2]*sqrt1mr2^3)*den
  fac  = 2r^2*onemr2*den
  facd = -2r*den
  t.dfdrbc[1] = t.c_n[1]*facd + t.c_n[2]*facd*sqrt1mr2
  for i=2:t.n
    flux -= t.c_n[i+1]*fac
    t.dfdrbc[1] += t.c_n[i+1]*facd*(2*onemr2-i*r^2)
    t.dfdrbc[i+3] -= fac
    fac *= sqrt1mr2
    facd *= sqrt1mr2
  end
  #  dfdrbc[2]=0 since the derivative with respect to b is zero.
  t.dfdrbc[3] = (onemr2-flux)*den
  t.dfdrbc[4] = 2/3*(sqrt1mr2^3-flux)*den
  # Also need to compute derivatives [ ]
  return flux
else
# Next, compute k^2 = m:
  onembmr2=(r-b+1)*(1-r+b); fourbr = 4b*r
  k2 = onembmr2/fourbr; 
  if k2 > 0
    k = sqrt(k2)
  else
    println("negative k2: ",k2," r: ",r," b: ",b)
  end
  dkdr = (b^2-r^2-1)/(8*k*b*r^2)
  dkdb = (r^2-b^2-1)/(8*k*b^2*r)
  if k2 > 1
    if k2 > 2.0
      kc = sqrt(1.-inv(k2))
    else
      kc2 = (1-r-b)*(1+b+r)/(1+r-b)/(1-r+b)
      kc = sqrt(kc2)
    end
  else
    if k2 > 0.5
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
    else
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
    end
  end
end

# Compute the highest value of v in J_v or I_v that we need:
# Compute sn[1] and its derivatives:
if b <= 1-r  # k^2 > 1
  lam = pi*r^2
  t.sn[1] = pi-lam
  t.dsndr[1] = -2*pi*r
  t.dsndb[1] = 0.
  kap0 = convert(T,pi); kck = zero(T)
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(sqarea_triangle(one(r),b,r)) 
  # Angle of section for occultor:
  kap0  = atan2(kite_area2,(r-1)*(r+1)+b^2)
  # Angle of section for source:
  pimkap1 = atan2(kite_area2,(r-1)*(r+1)-b^2)
  # Flux of visible uniform disk:
  t.sn[1] = pimkap1 - r^2*kap0 + .5*kite_area2
  t.dsndr[1]= -2*r*kap0
  t.dsndb[1]= kite_area2/b
  kck = kite_area2/(4*b*r)
end
# Compute sn[2] and its derivatives:
t.sn[2],Eofk,Em1mKdm = s2!(r,b,t.s2_grad)
t.dsndr[2] = t.s2_grad[1]
t.dsndb[2] = t.s2_grad[2]

# Compute the J_v and I_v functions:
if k2 > 0
  if (k2 < 0.5 || k2 > 2.0) # && v_max > 3
# This computes I_v,J_v for the largest v, and then works down to smaller values:
    dIJv_lower_dk!(k2,kck,kc,kap0,Eofk,Em1mKdm,t)
  else
# This computes I_0,J_0,J_1, and then works upward to larger v:
    dIJv_raise_dk!(k2,kck,kc,kap0,Eofk,Em1mKdm,t)
  end
end

# Next, loop over the Green's function components:
for n=2:t.n
  pofgn = zero(T)
  dpdr = zero(T)
  dpdb = zero(T)
  dpdk = zero(T)
  if iseven(n)
# For even values of n, sum over I_v:
    n0 = convert(Int64,n/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*t.Iv[n0+1]+2b*t.Iv[n0+2])
    dpdr = coeff*t.Iv[n0+1]
    dpdb = coeff*(-t.Iv[n0+1]+2*t.Iv[n0+2])
    dpdr += (n0+1)*pofgn/r
    dpdb += n0*pofgn/b
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
      term =  coeff*((r-b)*t.Iv[n0-i+1]+2b*t.Iv[n0-i+2])
      pofgn += term
      dpdr += coeff*t.Iv[n0-i+1]
      dpdb += coeff*(-t.Iv[n0-i+1]+2*t.Iv[n0-i+2])
      dpdr += term*(i*2*(b-r)/onembmr2+(n0+1-i)/r)
      dpdb += term*(i*2*(r-b)/onembmr2+(n0-i)/b)
    end
    pofgn *= 2r
    dpdr *= 2r
    dpdb *= 2r
  else
# Now do the same for odd N_c in sum over J_v:
    n0 = convert(Int64,(n-3)/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*t.Jv[n0+1]+2b*t.Jv[n0+2])
    dpdr = coeff*t.Jv[n0+1]
    dpdb = coeff*(-t.Jv[n0+1]+2*t.Jv[n0+2])
    dpdr  += pofgn*(3*(b-r)*r+(n0+1)*onembmr2)/(r*onembmr2)
    dpdb  += pofgn*(3*(r-b)*b+n0*onembmr2)/(b*onembmr2)
    dpdk = coeff*((r-b)*t.dJvdk[n0+1]+2b*t.dJvdk[n0+2])
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
      term = coeff*((r-b)*t.Jv[n0-i+1]+2b*t.Jv[n0-i+2])
      pofgn += term
      dpdr  +=  coeff*t.Jv[n0-i+1]
      dpdb  +=  coeff*(-t.Jv[n0-i+1]+2*t.Jv[n0-i+2])
      dpdr  += term*((i*2+3)*(b-r)*r+(n0+1-i)*onembmr2)/(r*onembmr2)
      dpdb  += term*((i*2+3)*(r-b)*b+(n0-i)*onembmr2)/(onembmr2*b)
      dpdk  += coeff*((r-b)*t.dJvdk[n0-i+1]+2b*t.dJvdk[n0-i+2])
    end
    pofgn *= 2r*onembmr2^1.5
    dpdr  *= 2r*onembmr2^1.5
    dpdb  *= 2r*onembmr2^1.5
    dpdk  *= 2r*onembmr2^1.5
  end
# Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
# boundary for n > 0.
# Compute sn[n]:
  t.sn[n+1] = -pofgn
  t.dsndr[n+1] = -(dpdr+dpdk*dkdr)
  t.dsndb[n+1] = -(dpdb+dpdk*dkdb)
end
#if typeof(r) == Float64
#  println("r: ",r," b: ",b," s2 error: ",convert(Float64,s2(big(r),big(b)))-sn[2])
#end
# That's it!
# Compute derivatives with respect to the coefficients:
den = inv(pi*(t.c_n[1]+2*t.c_n[2]/3))
flux = zero(T)
t.dfdrbc[1]=zero(T)  # Derivative with respect to r
t.dfdrbc[2]=zero(T)  # Derivative with respect to b
for i=0:n
  # derivatives with respect to the coefficients:
  t.dfdrbc[i+3]= t.sn[i+1]*den
  # total flux:
  flux += t.c_n[i+1]*t.dfdrbc[i+3]
  # derivatives with respect to r and b:
  t.dfdrbc[1] += t.c_n[i+1]*t.dsndr[i+1]*den
  t.dfdrbc[2] += t.c_n[i+1]*t.dsndb[i+1]*den
end
# Include derivatives with respect to first two c_n parameters:
t.dfdrbc[3] -= flux*den*pi
t.dfdrbc[4] -= flux*den*2pi/3
return flux
end
