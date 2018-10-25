# Computes K_v(k) and L_v(k) vectors from Agol & Luger (2018), along
# the derivatives with respect to k using recursion.
include("transit_structure.jl")

function dKLv_raise_dk!(t::Transit_Struct{T})  where {T <: Real}
kc = t.kc; k=t.k; k2 = t.k2; kc2 = t.kc2; kck=t.kck; kap=t.kap; Eofk=t.Eofk; Em1mKdm=t.Em1mKdm
# Compute K_v, L_v for 0 <= v <= v_max = l_max+2
# Define k:
# Iterate upwards in v:
v = t.v_max
# Compute I_v via upward iteration on v:
if k2 < 1
# First, compute value for v=0:
  t.Kv[1] = kap
# Try something else:
# Next, iterate upwards in v:
  f0 = kck
# Loop over v, computing K_v and L_v from higher v:
  @inbounds for v=1:t.v_max
    t.Kv[v+1]=((v-0.5)*t.Kv[v]+f0)/v
    f0 *= kc2
  end
  if t.grad
    # Now compute compute derivatives:
    t.dKvdk[1] = 2/kc  # TBD
    @inbounds for v=1:t.v_max
      t.dKvdk[v+1] = k2*t.dKvdk[v] # TB
    end
  end
else # k^2 >= 1
  # Compute v=0
  t.Kv[1] = pi
  @inbounds for v=1:t.v_max
    t.Kv[v+1]=t.Kv[v]*(1-1/(2v))
  end
  if t.grad
    # Derivatives of I_v are zero:
    fill!(t.dKvdk,zero(T))
  end
end
# Need to compute J_v for v=0, 1:
v= 0
if k2 < 1
  # Use elliptic integrals (already computed for s_2/S_1):
  if k2 > 0
    t.Lv[v+1]=2/(3k2*k)*(k2*(3k2-2)*Em1mKdm+k2*Eofk)
    t.Lv[v+2]= 2/(15k2*k)*(k2*(3k2+1)*Eofk-2k2*(1-3k2)*Em1mKdm)
    if t.grad
      t.dLvdk[v+1] = 2/k2*(-Eofk+2*Em1mKdm) # TBD
      t.dLvdk[v+2] = -3*t.Lv[v+2]/k+k2*t.dLvdk[v+1] # TBD
    end
  else
    t.Lv[v+1]= 0.0
    t.Lv[v+2]= 0.0
    if t.grad
      t.dLvdk[v+1] = 0.0
      t.dLvdk[v+2] = 0.0
    end
  end
else # k^2 >=1
  k2inv = inv(k2)
  t.Lv[v+1]=2/3*((3-2*k2inv)*Eofk+k2inv*Em1mKdm)
  t.Lv[v+2]=2*((3+k2inv)*Em1mKdm+(6-2*k2inv)*Eofk)/15
  if t.grad
    t.dLvdk[v+1] = 2/(k2*k)*(2*Eofk-Em1mKdm) # TBD
    t.dLvdk[v+2] = -3*t.Lv[v+2]/k+k2*t.dLvdk[v+1] # TBD
  end
end
@inbounds for v=2:t.v_max
#  t.Lv[v+1] = (2*(2v+(1-v)*k2)*t.Lv[v]+(k2-1)*(2v-3)*t.Lv[v-1])/(2v+3)
  t.Lv[v+1] = k2*(2*(2v*inv(k2)+(1-v))*t.Lv[v]+(1-inv(k2))*(2v-3)*t.Lv[v-1])/(2v+3)
  if t.grad
    t.dLvdk[v+1] = -3*t.Lv[v+1]/k+k2*t.dLvdk[v] # TBD
  end
end
return
end
