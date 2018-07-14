include("compute_c_n.jl")

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

#include("sn.jl")
include("s2.jl")
include("IJv_derivative.jl")
#include("area_triangle.jl")

#function compute_c_n(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
## Transform the u_n coefficients to c_n, which are coefficients
## of the basis in which the P(G_n) functions are computed.
#n = length(u_n)
#c_n = zeros(typeof(r),n+3)
#a_n = zeros(typeof(r),n+1)
## Pre-allocated memory - just need to zero terms:
#a_n[1] = one(r)  # Add in the first constant coefficient term
#for i=1:n
#  # Compute the contribution to a_n*\mu^n
#  bcoeff = one(T)
#  for j=0:i
##    a_n[j+1] -= u_n[i]*binomial(i,j)*(-1)^j
#    a_n[j+1] -= u_n[i]*bcoeff*(-1)^j
#    bcoeff *= (i-j)/(j+1)
##    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
#  end
#end
## Now, compute the c_n coefficients:
#for j=n:-1:2
#  c_n[j+1] = a_n[j+1]/(j+2)+c_n[j+3]
#end
#c_n[2] = a_n[2]+3*c_n[4]
#c_n[1] = a_n[1]+2*c_n[3]
#return c_n[1:n+1]
#end

function transit_poly_c(r::T,b::T,c_n::Array{T,1}) where {T <: Real}
# Number of limb-darkening components to include (beyond 0 and 1):
N_c = length(c_n)-1
# We are parameterizing these with the function:
# g_n = c_n [(n+2) z^n - n z^{n-2}] for n >= 2
# while g_{0} = c_0 z^0 (uniform source) and g_{1} = c_1 z^1 (linear limb-darkening)
# which gives a Green's integral of:
# P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
# for which we have a solution in terms of I_v (for even n) and J_v (for odd n).

# Set up a vector for storing results of P(G_n)-Q(G_n); note that
# this is a different vector than the Starry case:
#if ~prealloc
sn = zeros(typeof(r),N_c+1)
#end

# Check for different cases:
if b >= 1+r
  # unobscured - return one:
  return one(r)
end
if r >= 1+b
  # full obscuration - return zero:
  return zero(r)
end
if b == 0.0
  # Annular eclipse - integrate around the full boundary of both bodies:
  flux = zero(r); sqrt1mr2 = sqrt(1-r^2)
  flux = (c_n[1]*(1-r^2)+2/3*c_n[2]*sqrt1mr2^3)
  fac= 2r^2*(1-r^2)
  for i=2:N_c
    flux += -c_n[i+1]*fac
    fac *= sqrt1mr2
  end
  return flux/(c_n[1]+2*c_n[2]/3)
else
# Next, compute k^2 = m:
  onembmr2=(r+1-b)*(1-r+b); fourbr = 4b*r
#  onembpr2 = (1-r-b)*(1+b+r); onembmr2=(r-b+1)*(1-r+b); fourbr = 4b*r
#  k2 = onembpr2/fourbr+1
  k2 = onembmr2/fourbr
  if k2 > 1
    if k2 > 2.0
      kc = sqrt(1.-inv(k2))
    else
#      kc2 = (1-(b+r)^2)/(1-(b-r)^2)
      kc2 = (1-r-b)*(1+b+r)/(1-b+r)/(1-r+b)
      kc = sqrt(kc2)
    end
  else
    if k2 > 0.5
#      kc2 = ((b+r)^2-1)/(4*b*r)
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
    else
      kc = sqrt(1.-k2)
    end
  end
end

# Compute the highest value of v in J_v or I_v that we need:
if iseven(N_c)
  v_max = round(Int64,N_c/2)+2
else
  v_max = round(Int64,(N_c-1)/2)+2
end
# Compute sn[1] and sn[2]:
# Uniform disk case:
if b <= 1-r
  lam = pi*r^2
  sn[1] = pi-lam
  kap0 = convert(typeof(r),pi); kck = zero(r)
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(sqarea_triangle([one(r),b,r]))
  # Angle of section for occultor:
  kap0  = atan2(kite_area2,(r-1)*(r+1)+b^2)
  # Angle of section for source:
  pimkap1 = atan2(kite_area2,(r-1)*(r+1)-b^2)
  # Flux of visible uniform disk:
  sn[1] = pimkap1 - r^2*kap0 + .5*kite_area2
  kck = kite_area2/(4*b*r)
end
sn[2] = s2(r,b)
#if typeof(r) == Float64
# Compute the J_v and I_v functions:
#if ~prealloc
Iv = zeros(typeof(k2),v_max+1); Jv = zeros(typeof(k2),v_max+1)
#end
if k2 > 0
  if (k2 < 0.5 || k2 > 2.0) # && v_max > 3
# This computes I_v,J_v for the largest v, and then works down to smaller values:
    IJv_lower!(v_max,k2,kck,kc,kap0,Iv,Jv)
  else
# This computes I_0,J_0,J_1, and then works upward to larger v:
    IJv_raise!(v_max,k2,kck,kc,kap0,Iv,Jv)
  end
end

#nphi = 1000; dphi=2pi/nphi; phigrid = linspace(.5*dphi,1-.5*dphi,nphi)
# Next, loop over the Green's function components:
for n=2:N_c
  pofgn = zero(r)
  if iseven(n)
# For even values of n, sum over I_v:
    n0 = convert(Int64,n/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*Iv[n0+1]+2b*Iv[n0+2])
# For even n, compute coefficients for the sum over I_v:
#    println("n0: ",n0," i: ",0," coeff: ",coeff)
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
#      println("n0: ",n0," i: ",i," coeff: ",coeff)
      pofgn += coeff*((r-b)*Iv[n0-i+1]+2b*Iv[n0-i+2])
    end
    pofgn *= 2r
  else
# Now do the same for odd N_c in sum over J_v:
    n0 = convert(Int64,(n-3)/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*Jv[n0+1]+2b*Jv[n0+2])
#    println("n0: ",n0," i: ",0," coeff: ",coeff)
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
#      println("n0: ",n0," i: ",i," coeff: ",coeff)
      pofgn += coeff*((r-b)*Jv[n0-i+1]+2b*Jv[n0-i+2])
    end
    pofgn *= 2r*onembmr2^1.5
  end
#  pofgn_num = sum(sqrt.(1-r^2-b^2-2*b*r*sin.(phigrid)).^n.*(r+b.*sin.(phigrid))*r*dphi)
#  println("n: ",n," P(G_n): ",pofgn," P(G_n),num: ",pofgn_num)
# Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
# boundary for n > 0.
# Compute sn[n]:
  #println("n: ",n," P(G_n): ",pofgn)
  sn[n+1] = -pofgn
end
#  println("r: ",r," b: ",b," s2 error: ",convert(Float64,s2(big(r),big(b)))-sn[2])
#end
# That's it!
#println("s_n: ",sn)
#println("c_n*s_n: ",c_n.*sn)
flux = sum(c_n.*sn)/(pi*(c_n[1]+2*c_n[2]/3))  # for c_2 and above, the flux is zero.
return flux
end
# That's it!

# Now, for the versions of these functions which include derivatives:
#function compute_c_n(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
## Transform the u_n coefficients to c_n, which are coefficients
## of the basis in which the P(G_n) functions are computed.
## Compute the derivatives of the flux with respect to the u coefficients.
#n = length(u_n)
## We define c_n with two extra elements which are zero:
#c_n = zeros(typeof(r),n+3)
#dfdrbc = zeros(typeof(r),n+3)
#a_n = zeros(typeof(r),n+1)
#dadu = zeros(typeof(r),n+1,n)
#dcdu = zeros(typeof(r),n+3,n)
#a_n[1] = one(r)  # Add in the first constant coefficient term
#for i=1:n
#  # Compute the contribution to a_n*\mu^n
#  bcoeff = one(T)
#  for j=0:i
##    a_n[j+1] -= u_n[i]*binomial(i,j)*(-1)^j
#    a_n[j+1] -= u_n[i]*bcoeff*(-1)^j
##    dadu[j+1,i] -= binomial(i,j)*(-1)^j
#    dadu[j+1,i] -= bcoeff*(-1)^j
#    bcoeff *= (i-j)/(j+1)
##    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
#  end
#end
## Now, compute the c_n coefficients and propagate derivatives:
#for j=n:-1:2
#  c_n[j+1] = a_n[j+1]/(j+2)+c_n[j+3]
#  for i=1:n
#    dcdu[j+1,i] = dadu[j+1,i]/(j+2) + dcdu[j+3,i]
#  end
#end
#c_n[2] = a_n[2]+3*c_n[4]
#for i=1:n
#  dcdu[2,i] = dadu[2,i] + 3*dcdu[4,i]
#end
#c_n[1] = a_n[1]+2*c_n[3]
#for i=1:n
#  dcdu[1,i] = dadu[1,i] + 2*dcdu[3,i]
#end
#return c_n[1:n+1];dcdu[1:n+1,n]
#end

function transit_poly(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
# Compute the c_n values from u_n:
c_n = compute_c_n(u_n)
# Pass c_n (without last two dummy values):
flux = transit_poly_c(r,b,c_n)
return flux
end

function transit_poly!(r::T,b::T,c_n::Array{T,1},dcdu::Array{T,2},dfdrbc::Array{T,1},dfdrbu::Array{T,1}) where {T <: Real}
# Pass c_n (without last two dummy values):
flux = transit_poly_c!(r,b,c_n,dfdrbc)
# Now, transform derivaties from c to u:
fill!(dfdrbu,zero(r))
dfdrbu[1] = dfdrbc[1]  # r derivative
dfdrbu[2] = dfdrbc[2]  # b derivative
# u_n derivatives:
for i=1:n, j=0:n
  dfdrbu[i+2] += dfdrbc[j+3]*dcdu[j+1,i]
end
#println("dcdu: ",dcdu," dfdrbc: ",dfdrbc," dfdrbu: ",dfdrbu)
return flux
end

function transit_poly_c!(r::T,b::T,c_n::Array{T,1},dfdrbc::Array{T,1}) where {T <: Real}
@assert((length(c_n)+2) == length(dfdrbc))
@assert(r > 0)
# Number of limb-darkening components to include (beyond 0 and 1):
N_c = length(c_n)-1
# We are parameterizing these with the function:
# g_n = c_n [(n+2) z^n - n z^{n-2}] for n >= 2
# while g_{0} = c_0 z^0 (uniform source) and g_{1} = c_1 z^1 (linear limb-darkening)
# which gives a Green's integral of:
# P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
# for which we have a solution in terms of I_v (for even n) and J_v (for odd n).
# Compute the derivative of the flux with respect to the different coefficients.

# Set up a vector for storing results of P(G_n)-Q(G_n); note that
# this is a different vector than the Starry case:
#if ~prealloc
#  sn = zeros(typeof(r),N_c+1)
#  dsndr = zeros(typeof(r),N_c+1)
#  dsndb = zeros(typeof(r),N_c+1)
#end
# Check for different cases:
if b >= 1+r || r ==  0.0
  # unobscured - return one:
  return one(r)
end
if r >= 1+b
  # full obscuration - return zero:
  return zero(r)
end
if b == 0.0
  # Annular eclipse - integrate around the full boundary of both bodies:
  flux = zero(r); onemr2 = 1-r^2; sqrt1mr2 = sqrt(onemr2)
  den = inv(c_n[1]+2*c_n[2]/3)
  flux = (c_n[1]*onemr2+2/3*c_n[2]*sqrt1mr2^3)*den
  fac  = 2r^2*onemr2*den
  facd = -2r*den
  dfdrbc[1] = c_n[1]*facd + c_n[2]*facd*sqrt1mr2
  for i=2:N_c
    flux -= c_n[i+1]*fac
    dfdrbc[1] += c_n[i+1]*facd*(2*onemr2-i*r^2)
    dfdrbc[i+3] -= fac
    fac *= sqrt1mr2
    facd *= sqrt1mr2
  end
  #  dfdrbc[2]=0 since the derivative with respect to b is zero.
  dfdrbc[3] = (onemr2-flux)*den
  dfdrbc[4] = 2/3*(sqrt1mr2^3-flux)*den
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
#  onembpr2 = (1-b-r)*(1+b+r); onembmr2=(r-b+1)*(1-r+b); fourbr = 4b*r
#  k2 = onembpr2/fourbr+1; k = sqrt(k2)
  dkdr = (b^2-r^2-1)/(8*k*b*r^2)
#  dkdr = ((b-r)*(b+r)-1)/(8*k*b*r^2)
  dkdb = (r^2-b^2-1)/(8*k*b^2*r)
#  dkdb = ((r-b)*(r+b)-1)/(8*k*b^2*r)
  if k2 > 1
    if k2 > 2.0
      kc = sqrt(1.-inv(k2))
    else
#      kc2 = (1-(b+r)^2)/(1-(b-r)^2)
      kc2 = (1-r-b)*(1+b+r)/(1+r-b)/(1-r+b)
      kc = sqrt(kc2)
    end
#    if typeof(kc) == Float64
#      kc_error = kc-convert(Float64,sqrt((big(1.0)-big(b)-big(r))*(big(1.0)+big(b)+big(r))/(big(1.0)-big(b)+big(r))/(big(1.0)-big(r)+big(b))))
#      println("kc: ",kc," error on kc: ",kc_error)
#    end
  else
    if k2 > 0.5
#      kc2 = ((b+r)^2-1)/(4*b*r)
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
    else
      kc2 = (r-1+b)*(b+r+1)/(4*b*r)
      kc = sqrt(kc2)
#      kc = sqrt(1.-k2)
    end
#    if typeof(kc) == Float64
#      kc_error = kc-convert(Float64,sqrt((big(b)+big(r)-big(1.0))*(big(1.0)+big(b)+big(r))/(big(4.0)*big(b)*big(r))))
#      println("kc: ",kc," error on kc: ",kc_error)
#    end
  end
end

# Compute the highest value of v in J_v or I_v that we need:
if iseven(N_c)
  v_max = round(Int64,N_c/2)+2
else
  v_max = round(Int64,(N_c-1)/2)+2
end
#println("v_max: ",v_max," N_c: ",N_c)
# Compute sn[1] and its derivatives:
if b <= 1-r  # k^2 > 1
  lam = pi*r^2
  sn[1] = pi-lam
  dsndr[1] = -2*pi*r
  dsndb[1] = 0.
  kap0 = convert(typeof(r),pi); kck = zero(r)
else
  # Twice area of kite-shaped region connecting centers of circles & intersection points:
  kite_area2 = sqrt(sqarea_triangle([one(r),b,r])) 
  # Angle of section for occultor:
  kap0  = atan2(kite_area2,(r-1)*(r+1)+b^2)
  # Angle of section for source:
  pimkap1 = atan2(kite_area2,(r-1)*(r+1)-b^2)
  # Flux of visible uniform disk:
  sn[1] = pimkap1 - r^2*kap0 + .5*kite_area2
  dsndr[1]= -2*r*kap0
  dsndb[1]= kite_area2/b
  kck = kite_area2/(4*b*r)
end
# Compute sn[2] and its derivatives:
#if ~prealloc
#  s2_grad = zeros(typeof(r),2)
#end
sn[2] = s2!(r,b,s2_grad)
dsndr[2] = s2_grad[1]
dsndb[2] = s2_grad[2]

# Compute the J_v and I_v functions:
#if ~prealloc
#  Iv = zeros(typeof(k2),v_max+1); Jv = zeros(typeof(k2),v_max+1)
## And their derivatives with respect to k:
#  dIvdk = zeros(typeof(k2),v_max+1); dJvdk = zeros(typeof(k2),v_max+1)
#end
if k2 > 0
  if (k2 < 0.5 || k2 > 2.0) # && v_max > 3
# This computes I_v,J_v for the largest v, and then works down to smaller values:
    dIJv_lower_dk!(v_max,k2,kck,kc,kap0,Iv,Jv,dIvdk,dJvdk)
  else
# This computes I_0,J_0,J_1, and then works upward to larger v:
    dIJv_raise_dk!(v_max,k2,kck,kc,kap0,Iv,Jv,dIvdk,dJvdk)
  end
end

#nphi = 1000; dphi=2pi/nphi; phigrid = linspace(.5*dphi,1-.5*dphi,nphi)
# Next, loop over the Green's function components:
for n=2:N_c
  pofgn = zero(r)
  dpdr = zero(r)
  dpdb = zero(r)
  dpdk = zero(r)
  if iseven(n)
# For even values of n, sum over I_v:
    n0 = convert(Int64,n/2)
    coeff = (-fourbr)^n0
#    dIv_fac = (1+(r-b)*(r+b))/(2r)
    # Compute i=0 term
    pofgn = coeff*((r-b)*Iv[n0+1]+2b*Iv[n0+2])
    dpdr = coeff*Iv[n0+1]
    dpdb = coeff*(-Iv[n0+1]+2*Iv[n0+2])
    dpdr += (n0+1)*pofgn/r
#    dpdr += coeff*((n0+2-(n0+1)*b/r)*Iv[n0+1]+2*(n0+1)*b/r*Iv[n0+2])
    dpdb += n0*pofgn/b
#    println("v: ",n0,"-Iv[v]: ",-Iv[n0+1]," 2Iv[v+1]: ",2Iv[n0+2]," diff: ",-Iv[n0+1]+2*Iv[n0+2])
#    dpdk = coeff*((r-b)*dIvdk[n0+1]+2b*dIvdk[n0+2])
#    dpdk = coeff*dIvdk[n0+1]*dIv_fac
# For even n, compute coefficients for the sum over I_v:
#    println("n0: ",n0," i: ",0," coeff: ",coeff)
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
#      println("n0: ",n0," i: ",i," coeff: ",coeff)
      term =  coeff*((r-b)*Iv[n0-i+1]+2b*Iv[n0-i+2])
      pofgn += term
      dpdr += coeff*Iv[n0-i+1]
      dpdb += coeff*(-Iv[n0-i+1]+2*Iv[n0-i+2])
      dpdr += term*(i*2*(b-r)/onembmr2+(n0+1-i)/r)
#      dpdr += coeff*((n0+1)+i*((b-r)*(r+b)-1)/onembmr2)*((1-b/r)*Iv[n0-i+1]+2b/r*Iv[n0-i+2])
      dpdb += term*(i*2*(r-b)/onembmr2+(n0-i)/b)
#      println("v: ",n0-i,"-Iv[v]: ",-Iv[n0-i+1]," 2Iv[v+1]: ",2Iv[n0-i+2]," diff: ",-Iv[n0-i+1]+2*Iv[n0-i+2])
#      dpdk += coeff*((r-b)*dIvdk[n0-i+1]+2b*dIvdk[n0-i+2])
#      dpdk += coeff*dIvdk[n0-i+1]*dIv_fac
#      dpdk += coeff*2*i/k*((r-b)*Iv[n0-i+1]+2b*Iv[n0-i+2])
    end
    pofgn *= 2r
    dpdr *= 2r
#    dpdr += (n0+1)*pofgn/r
#    dpdr += pofgn/r
    dpdb *= 2r
#    dpdb += n0*pofgn/b
#    dpdk *= 2r
  else
# Now do the same for odd N_c in sum over J_v:
    n0 = convert(Int64,(n-3)/2)
    coeff = (-fourbr)^n0
    # Compute i=0 term
    pofgn = coeff*((r-b)*Jv[n0+1]+2b*Jv[n0+2])
    dpdr = coeff*Jv[n0+1]
    dpdb = coeff*(-Jv[n0+1]+2*Jv[n0+2])
#    dpdr += n0*pofgn/r
#    dpdb += n0*pofgn/b
#    dpdr  += pofgn*(3*(b-r)/onembmr2+(n0+1)/r)
    dpdr  += pofgn*(3*(b-r)*r+(n0+1)*onembmr2)/(r*onembmr2)
#    dpdb  += pofgn*(3*(r-b)/onembmr2+n0/b)
    dpdb  += pofgn*(3*(r-b)*b+n0*onembmr2)/(b*onembmr2)
    dpdk = coeff*((r-b)*dJvdk[n0+1]+2b*dJvdk[n0+2])
#    println("n0: ",n0," i: ",0," coeff: ",coeff)
# For even n, compute coefficients for the sum over I_v:
    for i=1:n0
      coeff *= -(n0-i+1)/i*k2
#      println("n0: ",n0," i: ",i," coeff: ",coeff)
      term = coeff*((r-b)*Jv[n0-i+1]+2b*Jv[n0-i+2])
      pofgn += term
      dpdr  +=  coeff*Jv[n0-i+1]
      dpdb  +=  coeff*(-Jv[n0-i+1]+2*Jv[n0-i+2])
#      dpdr  += term*((i*2+3)*(b-r)/onembmr2+(n0+1-i)/r)
      dpdr  += term*((i*2+3)*(b-r)*r+(n0+1-i)*onembmr2)/(r*onembmr2)
#      dpdb  += term*((i*2+3)*(r-b)/onembmr2+(n0-i)/b)
      dpdb  += term*((i*2+3)*(r-b)*b+(n0-i)*onembmr2)/(onembmr2*b)
      dpdk  += coeff*((r-b)*dJvdk[n0-i+1]+2b*dJvdk[n0-i+2])
#      dpdk  += coeff*2*i/k*((r-b)*Jv[n0-i+1]+2b*Jv[n0-i+2])
    end
    pofgn *= 2r*onembmr2^1.5
    dpdr  *= 2r*onembmr2^1.5
    dpdb  *= 2r*onembmr2^1.5
#    dpdr += ((n0+1)/r+3*(b-r)/onembmr2)*pofgn
#    dpdr  += (1/r+3*(b-r)/onembmr2)*pofgn
#    dpdb += (n0/b-3*(b-r)/onembmr2)*pofgn
#    dpdb  += 3*(r-b)/onembmr2*pofgn
#    dpdb += (n*onembmr2-3*(b-r)*(r+b)-3)/(2*b*onembmr2)*pofgn
    dpdk  *= 2r*onembmr2^1.5
  end
#  pofgn_num = sum(sqrt.(1-r^2-b^2-2*b*r*sin.(phigrid)).^n.*(r+b.*sin.(phigrid))*r*dphi)
#  println("n: ",n," P(G_n): ",pofgn," P(G_n),num: ",pofgn_num)
# Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
# boundary for n > 0.
# Compute sn[n]:
  sn[n+1] = -pofgn
#  println("n: ",n," sn: ",sn[n+1]," dpdr: ",dpdr," dpdk*dkdr: ",dpdk*dkdr," dpdb: ",dpdb," dpdk*dkdb: ",dpdk*dkdb)
  dsndr[n+1] = -(dpdr+dpdk*dkdr)
  dsndb[n+1] = -(dpdb+dpdk*dkdb)
end
#if typeof(r) == Float64
#  println("r: ",r," b: ",b," s2 error: ",convert(Float64,s2(big(r),big(b)))-sn[2])
#end
# That's it!
#println("s_n: ",sn)
#println("c_n*s_n: ",c_n.*sn)
# Compute derivatives with respect to the coefficients:
den = inv(pi*(c_n[1]+2*c_n[2]/3))
flux = zero(r)
dfdrbc[1]=zero(r)  # Derivative with respect to r
dfdrbc[2]=zero(r)  # Derivative with respect to b
for n=0:N_c
  # derivatives with respect to the coefficients:
  dfdrbc[n+3]= sn[n+1]*den
  # total flux:
  flux += c_n[n+1]*dfdrbc[n+3]
  # derivatives with respect to r and b:
  dfdrbc[1] += c_n[n+1]*dsndr[n+1]*den
  dfdrbc[2] += c_n[n+1]*dsndb[n+1]*den
end
# Include derivatives with respect to first two c_n parameters:
dfdrbc[3] -= flux*den*pi
dfdrbc[4] -= flux*den*2pi/3
#flux = sum(c_n.*sn)*den   # for c_2 and above, the flux integrated over the star is zero.
return flux
end
