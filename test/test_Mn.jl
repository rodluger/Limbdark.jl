2
# Tests the code for computing the derivatives
# of M_n with respect to k.
# include("../src/transit_poly_struct.jl")
# include("../src/Mn_compute.jl")
# include("../src/transit_structure.jl")
# include("../src/cel_bulirsch.jl")
# include("random.jl")

# using PyPlot

# Randomizer seed
#using Random
seed = Random.seed!(42)

using QuadGK

@testset "Mn_test" begin

function Mn_num(k2::T,m::Int64) where {T <: Real}
# Numerically integrates M_n(k^2)/(4br)^m (note: m can be a half-integer).
# See 11/10/2018 notes.
f(x) = sqrt(k2-sin(x)^2)^m
if k2 < 1.0
  kap2 = convert(T,asin(sqrt(big(k2))))
  Mn,error = quadgk(f,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Mn,error = quadgk(f,-pi2,pi2,rtol=1e-15)
end
return Mn
end

# Carries out a test of M_n for k^2 and n_max:
function test_Mn(r::T,b::T) where {T <: Real}
n_max = 30
#n_max = 10
# Set up transit structures for Float64 and BigFloat:
u = ones(n_max+2)
prec_frac = zeros(n_max+3)
prec_abs  = zeros(n_max+3)
# Initialize the transit structure to pass to routines:
t = transit_init(r,b,u,true)  # Float64
t_big = transit_init(big(r),big(b),big.(u),false) # BigFloat
# Call the transit routines, thereby computing Mn:
transit_poly!(t)
transit_poly!(t_big)

# Now compute with numerical integration, and compare:
reltol = 1e-6
abstol = 1e-15
# Compare with numerical integration:
for m=0:t.n_max
  Mnn = Mn_num(t.k2,m)*(2*t.sqbr)^m
  test1 = isapprox(Mnn,t.Mn[m+1],atol=abstol,rtol = reltol)
  if !test1
    println("m: ",m," k2: ",t.k2," Mn_num: ",Mnn," M_n: ",t.Mn[m+1]," diff: ",Mnn-t.Mn[m+1]," M_n[big]: ",convert(Float64,t_big.Mn[m+1]))
  end
  test2 = isapprox(t.Mn[m+1],convert(Float64,t_big.Mn[m+1]),atol=abstol,rtol=reltol)
  prec_frac[m+1] = t.Mn[m+1]/convert(Float64,t_big.Mn[m+1])-1.0
  prec_abs[m+1] = asinh(t.Mn[m+1])-asinh(convert(Float64,t_big.Mn[m+1]))
  if !test2
    println("m: ",m," k2: ",t.k2," Mn: ",t.Mn[m+1]," M_n_big: ",convert(Float64,t_big.Mn[m+1])," diff: ",t.Mn[m+1]-convert(Float64,t_big.Mn[m+1]))
  end
  @test test1
  @test test2
end
return prec_frac,prec_abs
end

epsilon = 1e-8
r=100.0; b0=[r-1+epsilon,r,r+1-epsilon]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    prec_frac,prec_abs = test_Mn(r,b)
#    semilogy(abs.(prec_frac))
    #semilogy(abs.(prec_abs),":")
    println("r: ",r," b: ",b," k2: ",k2," max error frac: ",maximum(abs.(prec_frac))," max error abs: ",maximum(abs.(prec_abs)))
  end
end

r=0.01; b0=[epsilon,r,1-r-epsilon,1-r+epsilon,1.0,r+1-epsilon]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    prec_frac,prec_abs= test_Mn(r,b)
#    semilogy(abs.(prec_frac))
    #semilogy(abs.(prec_abs),":")
    println("r: ",r," b: ",b," k2: ",k2," max error frac: ",maximum(abs.(prec_frac))," max error abs: ",maximum(abs.(prec_abs)))
  end
end

ntest = 10
for i=1:ntest
  r= 2.0 ; b = 0.25
  while r > 1 + b || r < b-1
    r = 2rand(seed); b=2rand(seed)
  end
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  prec_frac,prec_abs = test_Mn(r,b)
#  semilogy(abs.(prec_frac))
  #semilogy(abs.(prec_abs),":")
  println("r: ",r," b: ",b," k2: ",k2," max error frac: ",maximum(abs.(prec_frac))," max error abs: ",maximum(abs.(prec_abs)))
end
#axis([0,32,1e-17,1e-5])
end
