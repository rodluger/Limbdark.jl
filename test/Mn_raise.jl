include("../src/cel_bulirsch.jl")
include("../src/transit_structure.jl")
include("../src/transit_poly_struct.jl")
include("../src/Mn_compute.jl")
using QuadGK

function Mn_num(k2::T,m::T) where {T <: Real}
# Numerically integrates M_n(k^2)/(4br)^m (note: m can be a half-integer).
# See 11/10/2018 notes.
f(x) = (k2-sin(x)^2)^m
g(x) = (k2-sin(x)^2)^m*sin(x)^2
if k2 < 1.0
  kap2 = convert(T,asin(sqrt(big(k2))))
  Mn,error = quadgk(f,-kap2,kap2,rtol=1e-15)
  Nn,error = quadgk(g,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Mn,error = quadgk(f,-pi2,pi2,rtol=1e-15)
  Nn,error = quadgk(g,-pi2,pi2,rtol=1e-15)
end
return Mn,Nn
end

n = 30
t= transit_init(0.1,0.5,ones(n)/n,false)
nb = 7
rgrid = [0.01,0.1,1.0,10.0,100.0]
for r in rgrid
Mnr = zeros(t.n_max+1)
Nnr = zeros(t.n_max+1)
Mnl = zeros(t.n_max+1)
Nnl = zeros(t.n_max+1)
for i=1:nb
#r = 10.^(rand()*3-2); ; b=rand()*((1+r)-b0)+b0; fourbr = 4*b*r
  b0 = maximum([1e-8,r-1])
  b = b0+i/nb*(1+r-b0)
  t.r = r; t.b = b; t.fourbr = 4*b*r; t.sqbr = sqrt(b*r)
  t.sqarea = sqarea_triangle(1.0,r,b)
  t.onembmr2 = (1.0-(r-b)^2); t.onembpr2 = 1.0-(b+r)^2
  t.sqonembmr2 = sqrt(t.onembmr2)
  t.k2 = t.onembmr2/t.fourbr; 
  t.onemr2mb2 = 1-r^2-b^2
  if t.k2 < 1.0
    println("k2: ",t.k2)
    t.k2inv=inv(t.k2); t.kc = sqrt(1.0-t.k2); t.k = sqrt(t.k2); t.kap0 = 2*asin(t.k)
    t.Eofk = cel_bulirsch(t.k2,t.kc,1.0,1.0,1.0-t.k2)
    t.Em1mKdm = cel_bulirsch(t.k2,t.kc,1.0,1.0,0.0)
    t.kite_area2 = sqrt(abs(t.sqarea))
  else
    println("k2: ",t.k2)
    t.k2inv=inv(t.k2); t.kc = sqrt(1.0-t.k2inv); t.k = sqrt(t.k2); t.kap0 = pi
    t.Eofk = cel_bulirsch(t.k2inv,t.kc,1.0,1.0,1.0-t.k2inv)
    t.Em1mKdm = cel_bulirsch(t.k2inv,t.kc,1.0,1.0,0.0)
  end
  Mn_raise!(t)
  Mnr .= t.Mn
  Nn_raise!(t)
  Nnr .= t.Nn
  Mn_lower!(t)
  Nn_lower!(t)
  Mnl .= t.Mn
  Nnl .= t.Nn
  Mnn = zeros(t.n_max+1)
  # Now, compute the N_n functions: N_n = (4br)^{n/2} \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2\phi)^{n/2} \sin^2\phi d\phi
  for m=0:t.n_max
    Mn_n,Nn_num = Mn_num(t.k2,float(m/2))
    Mnn[m+1] = Mn_n*sqrt(t.fourbr)^m
    Nn_num *= sqrt(t.fourbr)^m
    println(Mnr[m+1]," ",Mnl[m+1]," ",Mnn[m+1]," ",Mnr[m+1]-Mnn[m+1]," ",Mnl[m+1]-Mnn[m+1]," ",Mnl[m+1]-Mnr[m+1])
    # Now use function which doesn't involve division by b:
    if t.k2 < 0.5
      println("n: ",m," N_n: ",Nnl[m+1]," N_n num: ",Nn_num)
    else
      println("n: ",m," N_n: ",Nnr[m+1]," N_n num: ",Nn_num)
    end
  end
  println("k2: ",t.k2," r: ",r," b: ",b," max(M_n-M_{n,num}): ",maximum(abs,Mnr-Mnn)," ",maximum(abs,Mnl-Mnn))
  read(STDIN,Char)
end
end
