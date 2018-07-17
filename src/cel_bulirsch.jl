# Computes the function cel(kc,p,a,b) from Bulirsch (1969).
function cel_bulirsch(k2::T,p::T,a::T,b::T) where {T <: Real}
@assert (k2 <= one(T))
ca = sqrt(eps(k2))
# Avoid undefined k2=1 case:
if k2 != 1.0
  kc = sqrt(1.0-k2)
else
  kc = eps(k2)
end
return cel_bulirsch(k2,kc,p,a,b)
end

# Version called with kc (this is to improve precision of computation):
function cel_bulirsch(k2::T,kc::T,p::T,a::T,b::T) where {T <: Real}
@assert (k2 <= one(T))
ca = sqrt(eps(k2))
# Avoid undefined k2=1 case:
if k2 == 1.0 || kc == 0.0
  kc = eps(k2)
end
# Initialize values:
ee = kc; m=1.0
if p > 0.0
  p = sqrt(p); pinv = inv(p); b *= pinv
else
  q=k2; g=1.0-p; f = g-k2
  q *= (b-a*p); ginv = inv(g); p=sqrt(f*ginv); a=(a-b)*ginv
  pinv = inv(p)
  b = -q*ginv^2*pinv+a*p
end
# Compute recursion:
f=a; a += b*pinv; g=ee*pinv; b += f*g; b +=b; p +=g; g=m; m += kc
iter = 0; itmax = 50
while abs(g-kc) > g*ca && iter < itmax
  kc=sqrt(ee)
  kc += kc 
  ee = kc*m
  f=a
  pinv = inv(p)
  a += b*pinv
  g=ee*pinv
  b += f*g
  b +=b
  p +=g
  g=m
  m += kc
  iter +=1
end
if iter == itmax
  println("k2 ",k2," kc ",kc," abs(g-kc) ",abs(g-kc)," g*ca ",g*ca," cel ",0.5*pi*(a*m+b)/(m*(m+p)))
end
return 0.5*pi*(a*m+b)/(m*(m+p))
end

# Vector version to attempt to improve speed for multiple elliptic integrals with same k_c:
function cel_bulirsch(k2::T,kc::T,p::T,a1::T,a2::T,a3::T,b1::T,b2::T,b3::T) where {T <: Real}
# This assumes first value of a and b uses p; the rest have p=1.
@assert (k2 <= one(T))
ca = sqrt(eps(k2))
# Avoid undefined k2=1 case:
if k2 == 1.0 || kc == 0.0
  kc = eps(k2)
end
# Initialize values:
ee = kc; m=1.0
if p > 0.0
  p = sqrt(p); pinv = inv(p); b1 *= pinv
else
  q=k2; g=1.0-p; f = g-k2
  q *= (b1-a1*p); ginv = inv(g); p=sqrt(f*ginv); a1=(a1-b1)*ginv
  pinv = inv(p)
  b1 = -q*ginv^2*pinv+a1*p
end
# Compute recursion:
f1=a1
# First compute the first integral with p:
a1 += b1*pinv; g=ee*pinv; b1 += f1*g; b1 *=2; p +=g; g=m;
# Next, compute the remainder with p = 1:
p1 = one(T); g1=ee
f2 = a2; f3 = a3
a2 += b2; b2 += f2*g1; b2 *=2
a3 += b3; b3 += f3*g1; b3 *=2
p1 +=g1
g1 = m
m += kc
iter = 0; itmax = 50
while (abs(g-kc) > g*ca || abs(g1-kc) > g1*ca) && iter < itmax
  kc = sqrt(ee)
  kc += kc
  ee = kc*m
  f1 = a1; f2=a2; f3=a3
  pinv = inv(p)
  pinv1 = inv(p1)
  a1 += b1*pinv
  a2 += b2*pinv1
  a3 += b3*pinv1
  g = ee*pinv
  g1= ee*pinv1
  b1 += f1*g
  b2 += f2*g1
  b3 += f3*g1
  b1 *= 2
  b2 *= 2
  b3 *= 2
  p  += g
  p1 += g1
  g  = m
  g1 = m
  m += kc
  iter +=1
end
if iter == itmax
  println("k2 ",k2," kc ",kc," abs(g-kc) ",abs(g-kc)," g*ca ",g*ca," cel ",0.5*pi*(a*m+b)/(m*(m+p)))
end
f1 = 0.5*pi*(a1*m+b1)/(m*(m+p))
f2 = 0.5*pi*(a2*m+b2)/(m*(m+p1))
f3 = 0.5*pi*(a3*m+b3)/(m*(m+p1))
return f1::T,f2::T,f3::T
end
