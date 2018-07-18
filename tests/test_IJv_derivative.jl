# Tests the code for computing the derivatives
# of I_v and J_v with respect to k.
include("../src/IJv_derivative_struct.jl")
#include("../src/dJv_seriesdk.jl")

function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}
# How to compute (twice) area squared of triangle with
# high precision (Goldberg 1991):
a,b,c=reverse(sort([a,b,c]))
area = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
return area
end

@testset "IJv_derivative" begin

function test_IJv_derivative(k2::Float64)
v_max = 20
Iv = zeros(Float64,v_max+1); Jv = zeros(Float64,v_max+1)
Iv_bigp = zeros(BigFloat,v_max+1); Jv_bigp = zeros(BigFloat,v_max+1)
Iv_bigm = zeros(BigFloat,v_max+1); Jv_bigm = zeros(BigFloat,v_max+1)
dIvdk = zeros(Float64,v_max+1); dJvdk = zeros(Float64,v_max+1)
dIvdk_big = zeros(BigFloat,v_max+1); dJvdk_big = zeros(BigFloat,v_max+1)
dIvdk2_num = zeros(Float64,v_max+1); dJvdk2_num = zeros(Float64,v_max+1)
# This computes I_v for the largest v, and then works down to smaller values:
dq = big(1e-18)
if k2 <= 1
  kc = convert(Float64,sqrt(1.-big(k2)))
  k = convert(Float64,sqrt(big(k2)))
  kck = convert(Float64,sqrt((1.-big(k2))*big(k2)))
  kap = convert(Float64,2*asin(sqrt(big(k2))))
else
  kc = convert(Float64,sqrt(1.-inv(big(k2))))
  k = convert(Float64,sqrt(big(k2)))
  kap = convert(typeof(r),pi); kck=zero(r)
end
if k2 > 0
  if k2 < 0.5 || k2 > 2.0
    dIJv_lower_dk!(v_max,k2,kck,kc,kap,Iv,Jv,dIvdk,dJvdk)
    if k2 <=1 
      k2_big = big(k2)+dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-big(k2)-dq)
      kap_big = 2*asin(k_big)
      dIJv_lower_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
      k2_big = big(k2)-dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-big(k2)+dq)
      kap_big = 2*asin(k_big)
      dIJv_lower_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
#      dIJv_lower_dk!(v_max,big(k2)-dq,sqrt((big(k2)-dq)*(1-big(k2)+dq)),sqrt(1-big(k2)+dq),Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
    else
      k2_big = big(k2)+dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-inv(big(k2)+dq))
      kap_big = big(0.0)
      dIJv_lower_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
#      dIJv_lower_dk!(v_max,big(k2)+dq,sqrt((big(k2)-dq)*(1-big(k2)+dq)),sqrt(1-inv(big(k2)+dq)),Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
      k2_big = big(k2)-dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-inv(big(k2)-dq))
      kap_big = big(0.0)
      dIJv_lower_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
#      dIJv_lower_dk!(v_max,big(k2)-dq,sqrt((big(k2)+dq)*(1-big(k2)-dq)),sqrt(1-inv(big(k2)-dq)),Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
    end
    for v=0:v_max
      dIvdk2_num[v+1] = convert(Float64,(Iv_bigp[v+1]-Iv_bigm[v+1])/(2dq))
      dJvdk2_num[v+1] = convert(Float64,(Jv_bigp[v+1]-Jv_bigm[v+1])/(2dq))
      test1 =  isapprox(dIvdk[v+1],dIvdk2_num[v+1]*2*k,atol = 1e-20)
      test2 = isapprox(dJvdk[v+1],dJvdk2_num[v+1]*2*k,atol = 1e-20) 
      @test isapprox(dIvdk[v+1],dIvdk2_num[v+1]*2*k,atol = 1e-20)
      @test isapprox(dJvdk[v+1],dJvdk2_num[v+1]*2*k,atol = 1e-20)
      if ~test1 || ~test2
        println("v: ",v," k2: ",k2," kc: ",kc," Iv: ",Iv[v+1]," Iv_big: ",convert(Float64,Iv_bigp[v+1])," Jv: ",Jv[v+1]," Jv_big: ",convert(Float64,Jv_bigp[v+1])," dIvdk: ",dIvdk[v+1]," dIvdk_num: ",dIvdk2_num[v+1]*2*k," diff: ",dIvdk[v+1]-dIvdk2_num[v+1]*2*k," dJvdk: ",dJvdk[v+1]," dJvdk_num: ",dJvdk2_num[v+1]*2*k," diff: ",dJvdk[v+1]-dJvdk2_num[v+1]*2*k)
      end
    end
  else
    dIJv_raise_dk!(v_max,k2,kck,kc,kap,Iv,Jv,dIvdk,dJvdk)
    if k2 <=1
      k2_big = big(k2)+dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-k2_big)
      kap_big = 2*asin(k_big)
      dIJv_raise_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
#      dIJv_raise_dk!(v_max,big(k2)+dq,sqrt((big(k2)+dq)*(1-big(k2)-dq)),sqrt(1-big(k2)-dq),Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
      k2_big = big(k2)-dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-k2_big)
      kap_big = 2*asin(k_big)
      dIJv_raise_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
#      dIJv_raise_dk!(v_max,big(k2)-dq,sqrt((big(k2)-dq)*(1-big(k2)+dq)),sqrt(1-big(k2)+dq),Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
    else
      k2_big = big(k2)+dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-inv(big(k2)+dq))
      kap_big = big(0.0)
      dIJv_raise_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
      k2_big = big(k2)-dq
      k_big = sqrt(k2_big)
      kc_big = sqrt(1-inv(big(k2)-dq))
      kap_big = big(0.0)
      dIJv_raise_dk!(v_max,k2_big,kc_big*k_big,kc_big,kap_big,Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
#      dIJv_raise_dk!(v_max,big(k2)+dq,sqrt((big(k2)+dq)*(1-big(k2)-dq)),sqrt(1-inv(big(k2)+dq)),Iv_bigp,Jv_bigp,dIvdk_big,dJvdk_big)
#      dIJv_raise_dk!(v_max,big(k2)-dq,sqrt((big(k2)-dq)*(1-inv(big(k2)-dq))),sqrt(1-inv(big(k2)-dq)),Iv_bigm,Jv_bigm,dIvdk_big,dJvdk_big)
    end
    for v=0:v_max
      dIvdk2_num[v+1] = convert(Float64,(Iv_bigp[v+1]-Iv_bigm[v+1])/(2dq))
      dJvdk2_num[v+1] = convert(Float64,(Jv_bigp[v+1]-Jv_bigm[v+1])/(2dq))
      test1 = isapprox(dIvdk[v+1],dIvdk2_num[v+1]*2*k,atol = 1e-20)
      test2 = isapprox(dJvdk[v+1],dJvdk2_num[v+1]*2*k,atol = 1e-20)
      @test isapprox(dIvdk[v+1],dIvdk2_num[v+1]*2*k,atol = 1e-20)
      @test isapprox(dJvdk[v+1],dJvdk2_num[v+1]*2*k,atol = 1e-20)
#      println("v: ",v," k2: ",k2," kc: ",kc," dIvdk2: ",dIvdk[v+1]/(2*k)," dIvdk2_num: ",dIvdk2_num[v+1]," dJvdk2: ",dJvdk[v+1]/(2k)," dJvdk2_num: ",dJvdk2_num[v+1])
#      println("v: ",v," k2: ",k2," kc: ",kc," dIvdk: ",dIvdk[v+1]," dIvdk_num: ",dIvdk2_num[v+1]*2*k," diff: ",dIvdk[v+1]-dIvdk2_num[v+1]*2*k," dJvdk: ",dJvdk[v+1]," dJvdk2_num: ",dJvdk2_num[v+1]*2*k," diff: ",dJvdk[v+1]-dJvdk2_num[v+1]*2*k)
      if ~test1 || ~test2
        println("v: ",v," k2: ",k2," kc: ",kc," Iv: ",Iv[v+1]," Iv_big: ",convert(Float64,Iv_bigp[v+1])," Jv: ",Jv[v+1]," Jv_big: ",convert(Float64,Jv_bigp[v+1])," dIvdk: ",dIvdk[v+1]," dIvdk_num: ",dIvdk2_num[v+1]*2*k," diff: ",dIvdk[v+1]-dIvdk2_num[v+1]*2*k," dJvdk: ",dJvdk[v+1]," dJvdk_num: ",dJvdk2_num[v+1]*2*k," diff: ",dJvdk[v+1]-dJvdk2_num[v+1]*2*k)
      end
    end
  end
end
return
end

r=100.0; b0=[r-1+1e-10,r,r+1-1e-10]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  test_IJv_derivative(k2)
end
r=0.01; b0=[1e-10,r,1-r-1e-10,1-r+1e-10,1.0,r+1-1e-10]
for i=1:length(b0)
  b=b0[i]
  k2 = (1-b+r)*(1+b-r)/(4*b*r)
  test_IJv_derivative(k2)
end
test_IJv_derivative(.5*rand())
test_IJv_derivative(.5+.5*rand())
test_IJv_derivative(inv(.5+.5*rand()))
test_IJv_derivative(inv(.5*rand()))
end
