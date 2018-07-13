# Optimizes the code for computing the derivatives
# of I_v and J_v with respect to k.
include("../src/IJv_derivative.jl")

function optimize_IJv_derivative(k2::Float64)
v_max = 20
Iv_lower = zeros(Float64,v_max+1); Jv_lower = zeros(Float64,v_max+1)
Iv_raise = zeros(Float64,v_max+1); Jv_raise = zeros(Float64,v_max+1)
dIvdk_lower = zeros(Float64,v_max+1); dJvdk_lower = zeros(Float64,v_max+1)
dIvdk_raise = zeros(Float64,v_max+1); dJvdk_raise = zeros(Float64,v_max+1)
Iv_big_lower = zeros(BigFloat,v_max+1); Jv_big_lower = zeros(BigFloat,v_max+1)
Iv_big_raise = zeros(BigFloat,v_max+1); Jv_big_raise = zeros(BigFloat,v_max+1)
dIvdk_big_lower = zeros(BigFloat,v_max+1); dJvdk_big_lower = zeros(BigFloat,v_max+1)
dIvdk_big_raise = zeros(BigFloat,v_max+1); dJvdk_big_raise = zeros(BigFloat,v_max+1)
# This computes I_v for the largest v, and then works down to smaller values:
if k2 <= 1
  kc = convert(Float64,sqrt(1.-big(k2)))
  k = convert(Float64,sqrt(big(k2)))
  kck = convert(Float64,sqrt((1.-big(k2))*big(k2)))
  kap = convert(Float64,2*asin(sqrt(big(k2))))
  k2_big = big(k2)
  k_big = sqrt(k2_big)
  kc_big = sqrt(1-big(k2))
  kck_big = k_big*kc_big
  kap_big = 2*asin(k_big)
else
  kc = convert(Float64,sqrt(1.-inv(big(k2))))
  k = convert(Float64,sqrt(big(k2)))
  kap = convert(typeof(k2),pi); kck=zero(k2)
  k2_big = big(k2)
  k_big = sqrt(k2_big)
  kc_big = sqrt(big(1.0)-inv(big(k2)))
  kck_big = big(0.0)
  kap_big = big(pi)
end
if k2 > 0
  dIJv_lower_dk!(v_max,k2,kck,kc,kap,Iv_lower,Jv_lower,dIvdk_lower,dJvdk_lower)
  dIJv_lower_dk!(v_max,k2_big,kck_big,kc_big,kap_big,Iv_big_lower,Jv_big_lower,dIvdk_big_lower,dJvdk_big_lower)
  dIJv_raise_dk!(v_max,k2,kck,kc,kap,Iv_raise,Jv_raise,dIvdk_raise,dJvdk_raise)
  dIJv_raise_dk!(v_max,k2_big,kck_big,kc_big,kap_big,Iv_big_raise,Jv_big_raise,dIvdk_big_raise,dJvdk_big_raise)
#  for v=0:v_max
#    println("v: ",v," k2: ",k2," kc: ",kc," Iv_lower: ",Iv_lower[v+1]," Iv_big_lower: ",convert(Float64,Iv_big_lower[v+1]),
#    " Iv_raise: ",Iv_raise[v+1]," Iv_big_raise: ",convert(Float64,Iv_big_raise[v+1]),
#    " Jv_lower: ",Jv_lower[v+1]," Jv_big_lower: ",convert(Float64,Jv_big_lower[v+1]),
#    " Jv_raise: ",Jv_raise[v+1]," Jv_big_raise: ",convert(Float64,Jv_big_raise[v+1]))
##      " dIvdk_lower: ",dIvdk[v+1]," dIvdk_num: ",dIvdk2_num[v+1]*2*k," diff: ",dIvdk[v+1]-dIvdk2_num[v+1]*2*k," dJvdk: ",dJvdk[v+1]," dJvdk_num: ",dJvdk2_num[v+1]*2*k," diff: ",dJvdk[v+1]-dJvdk2_num[v+1]*2*k)
#  end
end
return Iv_lower,Jv_lower,Iv_raise,Jv_raise,dIvdk_lower,dJvdk_lower,dIvdk_raise,dJvdk_raise,convert(Array{Float64,1},Iv_big_lower),convert(Array{Float64,1},Jv_big_lower),convert(Array{Float64,1},Iv_big_raise),convert(Array{Float64,1},Jv_big_raise),convert(Array{Float64,1},dIvdk_big_lower),convert(Array{Float64,1},dJvdk_big_lower),convert(Array{Float64,1},dIvdk_big_raise),convert(Array{Float64,1},dJvdk_big_raise)
end

nk2 = 100
v_max = 20
Iv_lower = zeros(Float64,nk2,v_max+1); Jv_lower = zeros(Float64,nk2,v_max+1)
Iv_raise = zeros(Float64,nk2,v_max+1); Jv_raise = zeros(Float64,nk2,v_max+1)
dIvdk_lower = zeros(Float64,nk2,v_max+1); dJvdk_lower = zeros(Float64,nk2,v_max+1)
dIvdk_raise = zeros(Float64,nk2,v_max+1); dJvdk_raise = zeros(Float64,nk2,v_max+1)
Iv_big_lower = zeros(Float64,nk2,v_max+1); Jv_big_lower = zeros(Float64,nk2,v_max+1)
Iv_big_raise = zeros(Float64,nk2,v_max+1); Jv_big_raise = zeros(Float64,nk2,v_max+1)
dIvdk_big_lower = zeros(Float64,nk2,v_max+1); dJvdk_big_lower = zeros(Float64,nk2,v_max+1)
dIvdk_big_raise = zeros(Float64,nk2,v_max+1); dJvdk_big_raise = zeros(Float64,nk2,v_max+1)
k2 = logspace(0,2,nk2)
for i=1:nk2
  Iv_lower[i,:],Jv_lower[i,:],Iv_raise[i,:],Jv_raise[i,:],dIvdk_lower[i,:],dJvdk_lower,dIvdk_raise[i,:],dJvdk_raise[i,:],Iv_big_lower[i,:],Jv_big_lower[i,:],Iv_big_raise[i,:],Jv_big_raise[i,:],dIvdk_big_lower[i,:],dJvdk_big_lower[i,:],dIvdk_big_raise[i,:],dJvdk_big_raise = optimize_IJv_derivative(k2[i])
end
