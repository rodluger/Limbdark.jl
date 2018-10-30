# Tests the code for computing the derivatives
# of I_v and J_v with respect to k.
include("../src/IJv_derivative_struct.jl")
include("../src/transit_structure.jl")
include("../src/cel_bulirsch.jl")
# include("random.jl")

# Randomizer seed
#using Random
seed = Random.seed!(42)

using QuadGK

@testset "IJv_derivative" begin

function Iv_num(k2::T,v::Int64) where {T <: Real}
# Numerically integrates I_v(k^2):
f(x) = sin(x)^(2v)
if k2 < 1.0
  kap2 = convert(Float64,asin(sqrt(big(k2))))
  Iv,error = quadgk(f,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Iv,error = quadgk(f,-pi2,pi2,rtol=1e-15)
end
return Iv
end

function Jv_num(k2::T,v::Int64) where {T <: Real}
# Numerically integrates L_v(k^2):
f(x) = sin(x)^(2v)*(1-sin(x)^2/k2)^1.5
if k2 < 1.0
  kap2 = convert(Float64,asin(sqrt(big(k2))))
  Jv,error = quadgk(f,-kap2,kap2,rtol=1e-15)
else
  pi2 = 0.5*pi
  Jv,error = quadgk(f,-pi2,pi2,rtol=1e-15)
end
return Jv
end

# Initialize all the variables needed for computing Iv, Jv:
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

# Initialize all the variables needed for computing Iv, Jv:
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
function test_IJv_derivative(k2::Float64)
@assert (k2 > 0)
v_max = 30
# Set up arrays to hold finite-difference derivatives:
dIvdk2_num = zeros(Float64,v_max+1); dJvdk2_num = zeros(Float64,v_max+1)
# Finite difference step:
dq = big(1e-18)
# Set up transit structures for Float64 and BigFloat.
# Since we only care about k, set r=b for simplicity:
r_big = inv(sqrt(2*big(k2))); b_big = r_big
r = convert(Float64,r_big); b=r; u = ones((v_max-2)*2); u_big = big.(u)
# Initialize the transit structure to pass to routines:
t = transit_init(r,b,u,true)  # Float64
t_big = transit_init(r_big,b_big,u_big,false) # BigFloat
t_bigm = transit_init(r_big,b_big,u_big,false) # BigFloat
t_bigp = transit_init(r_big,b_big,u_big,false) # BigFloat
t.k2 = k2
t_big.k2 = big(k2)
# Initialize variables needed for computing Iv/Jv:
initialize_big(t)
initialize(t_big)
# Compute the derivatives in Float64 precision:
if t.k2 < 0.5 || t.k2 > 2.0
  dIJv_lower_dk!(t)
  dIJv_lower_dk!(t_big)
else
  dIJv_raise_dk!(t)
  dIJv_raise_dk!(t_big)
end
reltol = 1e-6
abstol = 1e-15
# Compare with numerical integration:
for v=0:t.v_max
  Ivn = Iv_num(k2,v)
#  println("v: ",v," k2: ",k2," Iv_num: ",Ivn," Iv: ",t.Iv[v+1]," Iv_big: ",convert(Float64,t_big.Iv[v+1])," Eofk: ",t.Eofk," Em1mKdm: ",t.Em1mKdm)
  @test isapprox(Ivn,t.Iv[v+1],atol=abstol,rtol = reltol)
  @test isapprox(t.Iv[v+1],convert(Float64,t_big.Iv[v+1]),atol=abstol,rtol=reltol)
  Jvn = Jv_num(k2,v)
  println("v: ",v," k2: ",k2," Jv_num: ",Jvn," Jv: ",t.Jv[v+1]," Jv_big: ",convert(Float64,t_big.Jv[v+1]))
  @test isapprox(Jvn,t.Jv[v+1],atol=abstol,rtol=reltol)
  @test isapprox(t.Jv[v+1],convert(Float64,t_big.Jv[v+1]),atol=abstol,rtol=reltol)
end
return
# Now, compute finite differences in BigFloat precision:
t_bigp.k2 = big(k2)+dq
initialize(t_bigp)
if t.k2 < 0.5 || t.k2 > 2.0
  dIJv_lower_dk!(t_bigp)
else
  dIJv_raise_dk!(t_bigp)
end
t_bigm.k2 = big(k2)-dq
initialize(t_bigm)
if t.k2 < 0.5 || t.k2 > 2.0
  dIJv_lower_dk!(t_bigm)
else
  dIJv_raise_dk!(t_bigm)
end
for v=0:v_max
  dIvdk2_num[v+1] = convert(Float64,(t_bigp.Iv[v+1]-t_bigp.Iv[v+1])/(2dq))
  dJvdk2_num[v+1] = convert(Float64,(t_bigp.Jv[v+1]-t_bigp.Jv[v+1])/(2dq))
  test1 = isapprox(t.dIvdk[v+1],dIvdk2_num[v+1]*2*t.k,atol = 1e-20)
  test2 = isapprox(t.dJvdk[v+1],dJvdk2_num[v+1]*2*t.k,atol = 1e-20)
  @test isapprox(t.dIvdk[v+1],dIvdk2_num[v+1]*2*t.k,atol = 1e-20)
  @test isapprox(t.dJvdk[v+1],dJvdk2_num[v+1]*2*t.k,atol = 1e-20)
  if ~test1 || ~test2
    println("v: ",v," k2: ",t.k2," kc: ",t.kc," Iv: ",t.Iv[v+1]," Iv_big: ",convert(Float64,t_bigp.Iv[v+1]),
            " Jv: ",t.Jv[v+1]," Jv_big: ",convert(Float64,t_bigp.nv[v+1])," dIvdk: ",t.dIvdk[v+1],
            " dIvdk_num: ",dIvdk2_num[v+1]*2*t.k," diff: ",t.dIvdk[v+1]-dIvdk2_num[v+1]*2*t.k,
            " dJvdk: ",t.dJvdk[v+1]," dJvdk_num: ",dJvdk2_num[v+1]*2*t.k," diff: ",t.dJvdk[v+1]-dJvdk2_num[v+1]*2*t.k)
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
    test_IJv_derivative(k2)
  end
end
r=0.01; b0=[eps,r,1-r-eps,1-r+eps,1.0,r+1-eps]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  if k2 >= 0
    test_IJv_derivative(k2)
  end
end
ntest = 10
for i=1:ntest
  test_IJv_derivative(.5*rand(seed))
  test_IJv_derivative(.5+.5*rand(seed))
  test_IJv_derivative(inv(.5+.5*rand(seed)))
  test_IJv_derivative(inv(.5*rand(seed)))
end
end
