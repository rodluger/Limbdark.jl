# Tests the code for computing the derivatives
# of M_m with respect to k.
include("../src/transit_poly_struct.jl")
include("../src/Mm_compute.jl")
include("../src/transit_structure.jl")
include("../src/cel_bulirsch.jl")
# include("random.jl")

# Randomizer seed
#using Random
seed = Random.seed!(42)

using QuadGK

@testset "Mm_derivative" begin

function Mm_num(k2::T,m::Int64) where {T <: Real}
# Numerically integrates M_m(k^2)/(4br)^m (note: m can be a half-integer).
# See 11/10/2018 notes.
f(x) = sqrt(k2-sin(x)^2)^m
if k2 < 1.0
  kap2 = convert(T,asin(sqrt(big(k2))))
  Mm,error = quadgk(f,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Mm,error = quadgk(f,-pi2,pi2,rtol=1e-15)
end
return Mm
end

# Initialize all the variables needed for computing Mm:
function initialize_big(t::Transit_Struct{T}) where {T <: Real}
k2 = t.k2
t.sqonembmr2 = sqrt(1-(t.r-t.b)^2)
t.sqarea = sqarea_triangle(one(T),t.r,t.b)
t.sqbr = sqrt(t.r*t.b)
if k2 < 1
  t.kc2 = convert(T,1.0-big(k2))
  t.kc = convert(T,sqrt(1.0-big(k2)))
  t.k = convert(T,sqrt(big(k2)))
  t.kck = convert(T,sqrt((1.0-big(k2))*big(k2)))
  t.kap0 = convert(T,2*asin(sqrt(big(k2))))
  Piofk,Eofk_big,Em1mKdm_big = cel_bulirsch(big(k2),sqrt(1-big(k2)),zero(BigFloat), one(BigFloat),one(BigFloat),one(BigFloat), one(BigFloat),1-big(k2),zero(BigFloat))
  t.Eofk = convert(T,Eofk_big)
  t.Em1mKdm = convert(T,Em1mKdm_big)
else
  t.k2inv = convert(T,inv(big(k2)))
  t.kc2 = convert(T,1.0-inv(big(k2)))
  t.kc = convert(T,sqrt(1.0-inv(big(k2))))
  t.k = convert(T,sqrt(big(k2)))
  t.kap0 = convert(T,pi); kck=zero(T)
  Piofk,Eofk_big,Em1mKdm_big = cel_bulirsch(inv(big(k2)),sqrt(1-inv(big(k2))),zero(BigFloat), one(BigFloat),one(BigFloat),one(BigFloat), one(BigFloat),1-inv(big(k2)),zero(BigFloat))
  t.Eofk = convert(T,Eofk_big)
  t.Em1mKdm = convert(T,Em1mKdm_big)
end
return
end

# Initialize all the variables needed for computing Mm:
function initialize(t::Transit_Struct{T}) where {T <: Real}
k2 = t.k2
t.sqonembmr2 = sqrt(1-(t.r-t.b)^2)
t.sqarea = sqarea_triangle(one(T),t.r,t.b)
t.sqbr = sqrt(t.r*t.b)
if k2 < 1
  t.kc2 = convert(T,1.0-big(k2))
  t.kc = convert(T,sqrt(1.0-big(k2)))
  t.k = convert(T,sqrt(big(k2)))
  t.kck = convert(T,sqrt((1.0-big(k2))*big(k2)))
  t.kap0 = convert(T,2*asin(sqrt(big(k2))))
  Piofk,t.Eofk,t.Em1mKdm = cel_bulirsch(t.k2,t.kc,zero(T), one(T),one(T),one(T), one(T),t.kc2,zero(T))
else
  t.k2inv = convert(T,inv(big(k2)))
  t.kc2 = convert(T,1.0-t.k2inv)
  t.kc = convert(T,sqrt(1.0-inv(big(k2))))
  t.k = convert(T,sqrt(big(k2)))
  t.kap0 = convert(T,pi); kck=zero(T)
  Piofk,t.Eofk,t.Em1mKdm = cel_bulirsch(inv(t.k2),t.kc,zero(T), one(T),one(T),one(T), one(T),t.kc2,zero(T))
end
return
end

# Carries out a test of M_m for k^2 and m_max:
function test_Mm(k2::Float64)
@assert (k2 > 0)
m_max = 30
# Set up arrays to hold finite-difference derivatives:
# Finite difference step:
dq = big(1e-18)
# Set up transit structures for Float64 and BigFloat.
# Since we only care about k, set r=b for simplicity:
r_big = inv(sqrt(2*big(k2))); b_big = r_big
r = convert(Float64,r_big); b=r; u = ones(2*m_max+1); u_big = big.(u)
# Initialize the transit structure to pass to routines:
t = transit_init(r,b,u,true)  # Float64
t_big = transit_init(r_big,b_big,u_big,false) # BigFloat
t_bigm = transit_init(r_big,b_big,u_big,false) # BigFloat
t_bigp = transit_init(r_big,b_big,u_big,false) # BigFloat
t.k2 = k2
t_big.k2 = big(k2)
# Initialize variables needed for computing M_m:
initialize_big(t)
initialize(t_big)
# Compute the derivatives in Float64 precision:
if t.k2 < 0.5 || t.k2 > 2.0
  Mm_lower!(t)
  Mm_lower!(t_big)
else
  Mm_raise!(t)
  Mm_raise!(t_big)
end
reltol = 1e-6
abstol = 1e-15
# Compare with numerical integration:
for m=0:t.m_max
  Mmn = Mm_num(k2,m)
  @test isapprox(Mmn,t.Mm[m+1],atol=abstol,rtol = reltol)
  if !isapprox(Mmn,t.Mm[m+1],atol=abstol,rtol = reltol)
    println("m: ",m," Mm_num: ",Mmn," M_m: ",t.Mm[m+1]," diff: ",Mmn-t.Mm[m+1])
  end
  @test isapprox(t.Mm[m+1],convert(Float64,t_big.Mm[m+1]),atol=abstol,rtol=reltol)
  if !isapprox(t.Mm[m+1],convert(Float64,t_big.Mm[m+1]),atol=abstol,rtol=reltol)
    println("m: ",m," Mm: ",t.Mm[m+1]," M_m_big: ",convert(Float64,t_big.Mm[m+1])," diff: ",t.Mm[m+1]-convert(Float64,t_big.Mm[m+1]))
  end
end
return
end

eps = 1e-8
r=100.0; b0=[r-1+eps,r,r+1-eps]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    test_Mm(k2)
  end
end
r=0.01; b0=[eps,r,1-r-eps,1-r+eps,1.0,r+1-eps]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    test_Mm(k2)
  end
end
ntest = 10
for i=1:ntest
  test_Mm(.5*rand(seed))
  test_Mm(.5+.5*rand(seed))
  test_Mm(inv(.5+.5*rand(seed)))
  test_Mm(inv(.5*rand(seed)))
end
end
