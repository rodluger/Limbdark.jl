include("../src/cel_bulirsch.jl")
include("../src/transit_structure.jl")
include("../src/Mm_compute.jl")
using QuadGK

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
nb = 7
rgrid = [0.1,1.0,10.0]
for r in rgrid
Mmr = zeros(t.m_max+1)
Mml = zeros(t.m_max+1)
for i=1:nb
#r = 10.^(rand()*3-2); ; b=rand()*((1+r)-b0)+b0; fourbr = 4*b*r
  b0 = maximum([1e-8,r-1])
  b = b0+i/nb*(1+r-b0)
  t.r = r; t.b = b; t.fourbr = 4*b*r; t.sqbr = sqrt(b*r)
  t.onembmr2 = (1.0-(r-b)^2); t.onembpr2 = 1.0-(b+r)^2
  t.sqonembmr2 = sqrt(t.onembmr2)
  t.k2 = t.onembmr2/t.fourbr; 
  if t.k2 < 1.0
    println("k2: ",t.k2)
    t.k2inv=inv(t.k2); t.kc = sqrt(1.0-t.k2); t.k = sqrt(t.k2); t.kap0 = 2*asin(t.k)
  else
    t.k2inv=inv(t.k2); t.kc = sqrt(1.0-t.k2inv); t.k = sqrt(t.k2); t.kap0 = pi
  end
  Mm_raise!(t)
  Mmr .= t.Mm
  Mm_lower!(t)
  Mml .= t.Mm
  Mmn = zeros(t.m_max+1)
  for m=0:t.m_max
    Mmn[m+1] = Mm_num(t.k2,float(m/2))*sqrt(t.fourbr)^m
    println(Mmr[m+1]," ",Mml[m+1]," ",Mmn[m+1]," ",Mmr[m+1]-Mmn[m+1]," ",Mml[m+1]-Mmn[m+1]," ",Mml[m+1]-Mmr[m+1])
  end

  println("k2: ",t.k2," r: ",r," b: ",b," max(M_m-M_{m,num}): ",maximum(abs,Mmr-Mmn)," ",maximum(abs,Mml-Mmn))
  read(STDIN,Char)
end
end
