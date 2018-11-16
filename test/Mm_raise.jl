include("../src/cel_bulirsch.jl")
include("../src/transit_structure.jl")
using QuadGK

# Series evaluation of:
#    M_m/(4*b*r)^{m/2} = \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2{x})^{m/2} dx
# for the top four values of m = m_max-3 to m_max.
function Mm_series!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_m:
k2 = t.k2
# Computing leading coefficient (n=0):
if k2 < 1
  tol = eps(k2)
# Add leading term to M_m:
  Mm = t.Mm_coeff[1,:,1]
# Now, compute higher order terms until desired precision is reached:
  k2n = one(T)   # k^{2n}
  term = zeros(T,4)
  for n =1:t.nmax-1
    k2n *= k2
    term = k2n*t.Mm_coeff[1,:,n+1]
    Mm += term
    if maximum(abs.(term)) < tol
      break
    end
  end
  fac = t.sqonembmr2^(t.m_max-3)
  for j=4:-1:1
    Mm[j] *= fac
    fac *= t.sqonembmr2
  end
  return Mm*t.k
else # k^2 >= 1
  tol = eps(inv(k2))
  Mm = t.Mm_coeff[2,:,1]
  k2n = one(T); term = zeros(T,4)
  @inbounds for n = 1:t.nmax-1
    k2n *= t.k2inv
    term = k2n*t.Mm_coeff[2,:,n+1]
    Mm += term
    if maximum(abs.(term)) < tol
      break
    end
  end
  fac = t.sqonembmr2^(t.m_max-3)
  for j=4:-1:1
    Mm[j] *= fac
    fac *= t.sqonembmr2
  end
  return Mm
end
end

function Mm_raise(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; m_max=t.m_max; k2 = t.k2; kc = t.kc; k2inv = t.k2inv
Mm = t.Mm
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
Mm = zeros(T,m_max+1)
# Compute first four integrals:
if k2 < 1.0
  Mm[1] = t.kap0
  Mm[2] = 2*t.sqbr*2*k2*cel_bulirsch(k2,kc,1.0,1.0,0.0)
  Mm[3] = Mm[1]*(1.0-r^2-b^2)+sqrt((1-(r-b)^2)*((b+r)^2-1))
#  Mm[3] = Mm[1]*(1.0-r^2-b^2)+t.kite_area2
  mu = (4*k2-2); lam = (3*k2-2)*(k2-1)
  Mm[4] = (2*t.sqbr)^3*t.twothird*cel_bulirsch(k2,kc,1.0,lam+mu,lam+mu*(1.0-k2))
else
  Mm[1] = pi
  Mm[2] = 2*sqrt(1-(r-b)^2)*cel_bulirsch(k2inv,kc,1.0,1.0,1.0-k2inv)
  Mm[3] = pi*(1-r^2-b^2)
  mu = 2*(2-k2inv); lam = k2inv-1.0
  Mm[4] = t.sqonembmr2^3*t.twothird*cel_bulirsch(k2inv,kc,1.0,lam+mu,lam+mu*(1.0-k2inv))
end
for m=4:m_max
  Mm[m+1]=2*(1-1/m)*(1-r^2-b^2)*Mm[m-1]-(1-2/m)*t.onembmr2*t.onembpr2*Mm[m-3]
end
return Mm
end

function Mm_lower(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; m_max=t.m_max; k2 = t.k2
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
Mm = zeros(T,m_max+1)
# Compute top four integrals:
fac = (2*t.sqbr)^m_max
Mm[m_max+1] = fac*Mm_num(k2,0.5*m_max)
fac /= 2*t.sqbr
Mm[m_max]   = fac*Mm_num(k2,0.5*m_max-0.5)
fac /= 2*t.sqbr
Mm[m_max-1] = fac*Mm_num(k2,0.5*m_max-1.0)
fac /= 2*t.sqbr
Mm[m_max-2] = fac*Mm_num(k2,0.5*m_max-1.5)
# Compute series version:
t.Mm[m_max-2:m_max+1] = reverse(Mm_series!(t))
println("series vs. num: ",Mm[m_max-2:m_max+1]," ",t.Mm[m_max-2:m_max+1])
# Now iterate downwards:
for m=m_max-4:-1:0
  Mm[m+1]=((m+4)*Mm[m+5]-2*(m+3)*(1-r^2-b^2)*Mm[m+3])/((1-(b-r)^2)*((b+r)^2-1)*(m+2))
end
return Mm
end

function Mm_num(k2::T,m::T) where {T <: Real}
# Numerically integrates M_m(k^2)/(4br)^m (note: m can be a half-integer).
# See 11/10/2018 notes.
f(x) = (k2-sin(x)^2)^m
if k2 < 1.0
  kap2 = convert(T,asin(sqrt(big(k2))))
  Mm,error = quadgk(f,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Mm,error = quadgk(f,-pi2,pi2,rtol=1e-15)
end
return Mm
end

n = 30
t= transit_init(0.1,0.5,ones(n)/n,false)
nb = 3
rgrid = [0.1,1.0,10.0]
for r in rgrid
for i=1:nb
#r = 10.^(rand()*3-2); b0 = maximum([0,r-1]); b=rand()*((1+r)-b0)+b0; fourbr = 4*b*r
  b = abs(1-r)+i/nb*(1+r-abs(1-r))
  t.r = r; t.b = b; t.fourbr = 4*b*r; t.sqbr = sqrt(b*r)
  t.onembmr2 = (1.0-(r-b)^2); t.onembpr2 = 1.0-(b+r)^2
  t.sqonembmr2 = sqrt(t.onembmr2)
  t.k2 = t.onembmr2/t.fourbr; t.k2inv=inv(t.k2); t.kc = sqrt(1.0-t.k2); t.k = sqrt(t.k2); t.kap0 = 2*asin(t.k)
  Mmr = Mm_raise(t)
  Mml = Mm_lower(t)
  Mmn = zeros(t.m_max+1)
  for m=0:t.m_max
    Mmn[m+1] = Mm_num(t.k2,float(m/2))*sqrt(t.fourbr)^m
    println(Mmr[m+1]," ",Mml[m+1]," ",Mmn[m+1]," ",Mmr[m+1]-Mmn[m+1]," ",Mml[m+1]-Mmn[m+1])
  end

  println("k2: ",t.k2," r: ",r," b: ",b," max(M_m-M_{m,num}): ",maximum(abs,Mmr-Mmn)," ",maximum(abs,Mml-Mmn))
  read(STDIN,Char)
end
end
