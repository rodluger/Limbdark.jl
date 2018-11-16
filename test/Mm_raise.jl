include("../src/cel_bulirsch.jl")
using QuadGK


# Series evaluation of:
#    M_m/(4*b*r)^{m/2} = \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2{x})^{m/2} dx
function Mm_series(k2::T,m::Int64,Mm_coeff::Array{T,2}) where {T <: Real}
nmax

return Mm
end

# Coefficients for series expansion of M_{m_max},..., M_{m_max-3}
function Mm_coeff(m_max::Int64,Mmc::Array{T,2},dMmcdk::Array{T,2},nmax::Int64) where {T <: Real}
# First, compute m_max:
for i=0:nmax
  Mmc[4,i+1] = 
return Mmc
end

function Mm_raise(r::T,b::T,m_max::Int64) where {T <: Real}
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
Mm = zeros(T,m_max+1)
fourbr = 4*b*r
k2 = (1.0-(r-b)^2)/fourbr
# Compute first four integrals:
if k2 < 1.0
  kc = sqrt(1.0-k2)
  Mm[1] = 2*asin(sqrt(k2))
  Mm[2] = sqrt(fourbr)*2*k2*cel_bulirsch(k2,kc,1.0,1.0,0.0)
  Mm[3] = Mm[1]*(1.0-r^2-b^2)+sqrt((1-(r-b)^2)*((b+r)^2-1))
  mu = (4*k2-2); lam = (3*k2-2)*(k2-1)
  Mm[4] = sqrt(fourbr)^3*2/3*cel_bulirsch(k2,kc,1.0,lam+mu,lam+mu*(1.0-k2))
else
  kc = sqrt(1.0-1.0/k2); k2inv = inv(k2)
  Mm[1] = pi
  Mm[2] = 2*sqrt(1-(r-b)^2)*cel_bulirsch(k2inv,kc,1.0,1.0,1.0-k2inv)
  Mm[3] = pi*(1-r^2-b^2)
  mu = 2*(2-k2inv); lam = k2inv-1.0
  Mm[4] = sqrt(fourbr*k2)^3*2/3*cel_bulirsch(k2inv,kc,1.0,lam+mu,lam+mu*(1.0-k2inv))
end
for m=4:m_max
  Mm[m+1]=2*(1-1/m)*(1-r^2-b^2)*Mm[m-1]+(1-2/m)*(1-(b-r)^2)*((b+r)^2-1)*Mm[m-3]
end
return Mm
end

function Mm_lower(r::T,b::T,m_max::Int64) where {T <: Real}
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
Mm = zeros(T,m_max+1)
fourbr = 4*b*r
k2 = (1.0-(r-b)^2)/fourbr
# Compute top four integrals:
sqfourbr = sqrt(fourbr)
sqfourbrinv = inv(sqfourbr)
fac = sqfourbr^m_max
Mm[m_max+1] = fac*Mm_num(k2,0.5*m_max)
fac *= sqfourbrinv
Mm[m_max]   = fac*Mm_num(k2,0.5*m_max-0.5)
fac *= sqfourbrinv
Mm[m_max-1] = fac*Mm_num(k2,0.5*m_max-1.0)
fac *= sqfourbrinv
Mm[m_max-2] = fac*Mm_num(k2,0.5*m_max-1.5)
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

m_max = 30
nb = 20
for i=1:nb
#r = 10.^(rand()*3-2); b0 = maximum([0,r-1]); b=rand()*((1+r)-b0)+b0; fourbr = 4*b*r
r = 0.1; b = 1-r+i/nb*2r; fourbr = 4*b*r
k2 = (1.0-(r-b)^2)/fourbr
Mmr = Mm_raise(r,b,m_max)
Mml = Mm_lower(r,b,m_max)
Mmn = zeros(m_max+1)
for m=0:m_max
  Mmn[m+1] = Mm_num(k2,float(m/2))*sqrt(fourbr)^m
  println(Mmr[m+1]," ",Mml[m+1]," ",Mmn[m+1]," ",Mmr[m+1]-Mmn[m+1]," ",Mml[m+1]-Mmn[m+1])
end

println("k2: ",k2," r: ",r," b: ",b," max(M_m-M_{m,num}): ",maximum(abs,Mmr-Mmn)," ",maximum(abs,Mml-Mmn))
read(STDIN,Char)
end
