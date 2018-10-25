# Tests the code for computing the derivatives
# of K_v and L_v with respect to k.
include("../src/KLv_derivative_struct.jl")
include("../src/transit_structure.jl")
include("../src/cel_bulirsch.jl")

# Randomizer seed
using Random
Random.seed!(42)
using QuadGK

@testset "KLv_derivative" begin

function Kv_num(k2::T,v::Int64) where {T <: Real}
# Numerically integrates K_v(k^2):
f(x) = cos(x)^(2v)
if k2 < 1.0
  kap2 = convert(Float64,asin(sqrt(big(k2))))
  Kv,error = quadgk(f,-kap2,kap2,rtol=1e-13)
else
  pi2 = 0.5*pi
  Kv,error = quadgk(f,-pi2,pi2,rtol=1e-13)
end
return Kv
end

function Lv_num(k2::T,v::Int64) where {T <: Real}
# Numerically integrates L_v(k^2):
f(x) = cos(x)^(2v)*(1-sin(x)^2/k2)^1.5
if k2 < 1.0
  kap2 = convert(Float64,asin(sqrt(big(k2))))
  Lv,error = quadgk(f,-kap2,kap2,rtol=1e-8)
else
  pi2 = 0.5*pi
  Lv,error = quadgk(f,-pi2,pi2,rtol=1e-8)
end
return Lv
end

# Initialize all the variables needed for computing Kv, Lv:
function initialize_big(t::Transit_Struct{T}) where {T <: Real}
k2 = t.k2
if k2 < 1
  t.kc2 = convert(T,1.0-big(k2))
  t.kc = convert(T,sqrt(1.0-big(k2)))
  t.k = convert(T,sqrt(big(k2)))
  t.kck = convert(T,sqrt((1.0-big(k2))*big(k2)))
  t.kap = convert(T,2*asin(sqrt(big(k2))))
  Piofk,Eofk_big,Em1mKdm_big = cel_bulirsch(big(k2),sqrt(1-big(k2)),zero(BigFloat), one(BigFloat),one(BigFloat),one(BigFloat), one(BigFloat),1-big(k2),zero(BigFloat))
  t.Eofk = convert(T,Eofk_big)
  t.Em1mKdm = convert(T,Em1mKdm_big)
else
  t.kc2 = convert(T,1.0-inv(big(k2)))
  t.kc = convert(T,sqrt(1.0-inv(big(k2))))
  t.k = convert(T,sqrt(big(k2)))
  t.kap = convert(T,pi); kck=zero(T)
  Piofk,Eofk_big,Em1mKdm_big = cel_bulirsch(inv(big(k2)),sqrt(1-inv(big(k2))),zero(BigFloat), one(BigFloat),one(BigFloat),one(BigFloat), one(BigFloat),1-inv(big(k2)),zero(BigFloat))
  t.Eofk = convert(T,Eofk_big)
  t.Em1mKdm = convert(T,Em1mKdm_big)
end
return
end

# Initialize all the variables needed for computing Kv, Lv:
function initialize(t::Transit_Struct{T}) where {T <: Real}
k2 = t.k2
if k2 < 1
  t.kc2 = convert(T,1.0-big(k2))
  t.kc = convert(T,sqrt(1.0-big(k2)))
  t.k = convert(T,sqrt(big(k2)))
  t.kck = convert(T,sqrt((1.0-big(k2))*big(k2)))
  t.kap = convert(T,2*asin(sqrt(big(k2))))
  Piofk,t.Eofk,t.Em1mKdm = cel_bulirsch(t.k2,t.kc,zero(T), one(T),one(T),one(T), one(T),t.kc2,zero(T))
else
  t.kc2 = convert(T,1.0-inv(big(k2)))
  t.kc = convert(T,sqrt(1.0-inv(big(k2))))
  t.k = convert(T,sqrt(big(k2)))
  t.kap = convert(T,pi); kck=zero(T)
  Piofk,t.Eofk,t.Em1mKdm = cel_bulirsch(inv(t.k2),t.kc,zero(T), one(T),one(T),one(T), one(T),t.kc2,zero(T))
end
return
end

# Carries out a test of K_v and L_v for k^2 and v_max:
function test_KLv_derivative(k2::Float64)
@assert (k2 > 0)
v_max = 50
# Set up arrays to hold finite-difference derivatives:
dKvdk2_num = zeros(Float64,v_max+1); dLvdk2_num = zeros(Float64,v_max+1)
# Finite difference step:
dq = big(1e-18)
# Set up transit structures for Float64 and BigFloat.
# Since we only care about k, set r=b for simplicity:
r_big = inv(sqrt(2*big(k2))); b_big = r_big
r = convert(Float64,r_big); b=r; u = ones((v_max-2)*2); u_big = big.(u)
# Initialize the transit structure to pass to routines:
t = transit_init(r,b,u,true)  # Float64
t_bigm = transit_init(r_big,b_big,u_big,false) # BigFloat
t_big = transit_init(r_big,b_big,u_big,false) # BigFloat
t_bigp = transit_init(r_big,b_big,u_big,false) # BigFloat
t.k2 = k2
# Initialize variables needed for computing Kv/Lv:
initialize_big(t)
# Compute the derivatives in Float64 precision:
dKLv_raise_dk!(t)
dKLv_raise_dk!(t_big)
# Compare with numerical integration:
for v=0:t.v_max
  Kvn = Kv_num(k2,v)
  println("v: ",v," k2: ",k2," Kv_num: ",Kvn," Kv: ",convert(Float64,t_big.Kv[v+1]))
  @test isapprox(Kvn,convert(Float64,t_big.Kv[v+1]))
  Lvn = Lv_num(k2,v)
  println("v: ",v," k2: ",k2," Lv_num: ",Lvn," Lv: ",convert(Float64,t_big.Lv[v+1]))
  @test isapprox(Lvn,convert(Float64,t_big.Lv[v+1]))
end
return
# Now, compute finite differences in BigFloat precision:
t_bigp.k2 = big(k2)+dq
initialize(t_bigp)
dKLv_raise_dk!(t_bigp)
t_bigm.k2 = big(k2)-dq
initialize(t_bigm)
dKLv_raise_dk!(t_bigm)
for v=0:v_max
  dKvdk2_num[v+1] = convert(Float64,(t_bigp.Kv[v+1]-t_bigp.Kv[v+1])/(2dq))
  dLvdk2_num[v+1] = convert(Float64,(t_bigp.Lv[v+1]-t_bigp.Lv[v+1])/(2dq))
  test1 = isapprox(t.dKvdk[v+1],dKvdk2_num[v+1]*2*t.k,atol = 1e-20)
  test2 = isapprox(t.dLvdk[v+1],dLvdk2_num[v+1]*2*t.k,atol = 1e-20)
  @test isapprox(t.dKvdk[v+1],dKvdk2_num[v+1]*2*t.k,atol = 1e-20)
  @test isapprox(t.dLvdk[v+1],dLvdk2_num[v+1]*2*t.k,atol = 1e-20)
  if ~test1 || ~test2
    println("v: ",v," k2: ",t.k2," kc: ",t.kc," Kv: ",t.Kv[v+1]," Kv_big: ",convert(Float64,t_bigp.Kv[v+1]),
            " Lv: ",t.Lv[v+1]," Lv_big: ",convert(Float64,t_bigp.Lv[v+1])," dKvdk: ",t.dKvdk[v+1],
            " dKvdk_num: ",dKvdk2_num[v+1]*2*t.k," diff: ",t.dKvdk[v+1]-dKvdk2_num[v+1]*2*t.k,
            " dLvdk: ",t.dLvdk[v+1]," dLvdk_num: ",dLvdk2_num[v+1]*2*t.k," diff: ",t.dLvdk[v+1]-dLvdk2_num[v+1]*2*t.k)
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
    test_KLv_derivative(k2)
  end
end
r=0.01; b0=[eps,r,1-r-eps,1-r+eps,1.0,r+1-eps]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    test_KLv_derivative(k2)
  end
end
test_KLv_derivative(.5*rand())
test_KLv_derivative(.5+.5*rand())
test_KLv_derivative(inv(.5+.5*rand()))
test_KLv_derivative(inv(.5*rand()))
end
