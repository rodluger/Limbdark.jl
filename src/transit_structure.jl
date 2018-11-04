mutable struct Transit_Struct{T}
  # Structure to hold arrays and other quantities for computing transit:
  r       :: T           # radius ratio
  b       :: T           # impact parameter
  u_n     :: Array{T,1}  # limb-darkening coefficients
  n       :: Int64       # number of limb-darkening coefficients
  v_max   :: Int64       # maximum coefficient in I_v and J_v
  c_n     :: Array{T,1}  # Green's basis coefficients
  sn      :: Array{T,1}  # Green's basis terms
  Iv      ::  Array{T,1} # integral over sin(phi)^{2v}
  Jv      ::  Array{T,1} # integral over sin(phi)^{2v} (1-sin(phi)^2/k^2)^{3/2}
  Kv      ::  Array{T,1} # integral over cos(phi)^{2v}
  Lv      ::  Array{T,1} # integral over cos(phi)^{2v} (1-sin(phi)^2/k^2)^{3/2}
  grad    :: Bool        # true for gradient; false for no gradient
  dIvdk   ::  Array{T,1} # derivative of I_v with respect to k
  dJvdk   ::  Array{T,1} # derivative of J_v with respect to k
  dKvdk   ::  Array{T,1} # derivative of K_v with respect to k
  dLvdk   ::  Array{T,1} # derivative of L_v with respect to k
  s2_grad :: Array{T,1}  # gradient of s2 with respect to (r,b) - this is S_1 = s_2 from starry
  dsndr   ::  Array{T,1} # derivatives of Green's basis with respect to r
  dsndb   ::  Array{T,1} # derivatives of Green's basis with respect to b
  dcdu    :: Array{T,2}  # derivatives c_n with respect to u_n
  dfdrbc  :: Array{T,1}  # derivative of flux with respect to c_n
  dfdrbu  :: Array{T,1}  # derivative of flux with respect to u_n
  nmax    :: Int64       # maximum number of terms in series expansions of I_v and J_v
  Iv_coeff :: Array{T,1} # coefficients for series expansion of I_v
  Jv_coeff :: Array{T,3} # coefficients for series expansion of J_v
  dJvdk_coeff :: Array{T,3} # coefficients for series expansion of dJ_v/dk
  k2      :: T           # k^2 = (1-(r-b)^2)/(4br)
  k       :: T           # k = sqrt(k^2)
  kc      :: T           # k_c = sqrt(1-k^2) (unless k > 1, then it is sqrt(1-1/k^2))
  kck     :: T           # k_c * k
  kap0    :: T           # kappa = sin^{-1}(k) (=kappa_0 in M&A)
  pimkap1 :: T           # pi-kappa_1
  kite_area2 :: T        # 2*A_kite
  Eofk    :: T           # E(k^2) is complete elliptic integral of first kind
  Em1mKdm :: T           # (E(m)-(1-m)K(m))/m is complete elliptic integral with m=k^2
  onembmr2:: T           # 1-(b-r)^2
  onembpr2:: T           # 1-(b+r)^2
  onembmr2inv :: T       # 1/(1-(b-r)^2)
  sqonembmr2  :: T       # sqrt(1-(b-r)^2)
  fourbr  :: T           # 4*b*r
  sqbr    :: T           # sqrt(b*r)
  sqbrinv :: T           # 1/sqrt(b*r)
  fourbrinv  :: T           # 1/(4*b*r)
  k2inv   :: T           # 1/k^2
  kc2     :: T           # 1-k^2 (or 1-1/k^2 if k > 1)
  sqrt1mr2:: T           # sqrt(1-r^2)
  den     :: T           # 1/(c_1 + c_2*2/3)
  third   :: T           # 1/3
  twothird:: T           # 2/3
  sqr1mr  :: T           # sqrt(r*(1-r)) if r < 1
  bincoeff:: Array{T,2}  # binomial(n,i)
end

using SpecialFunctions
include("compute_c_n_struct.jl")
include("IJv_coeff.jl")

function transit_init(r::T,b::T,u_n::Array{T,1},grad::Bool) where {T <: Real}
# Initializs a transit structure.
n = length(u_n)
if iseven(n)
  v_max = round(Int64,n/2)+2
else
  v_max = round(Int64,(n-1)/2)+2
end
nmax = 50
trans = Transit_Struct{T}(r,b,u_n,n,v_max,
  zeros(T,n+1),    # c_n
  zeros(T,n+1),    # sn
  zeros(T,v_max+1),# Iv
  zeros(T,v_max+1),# Jv
  zeros(T,v_max+1),# Kv
  zeros(T,v_max+1),# Lv
  grad,            # grad
  zeros(T,v_max+1),# dIvdk
  zeros(T,v_max+1),# dJvdk
  zeros(T,v_max+1),# dKvdk
  zeros(T,v_max+1),# dLvdk
  zeros(T,2),      # s2_grad
  zeros(T,n+1),    # dsndr
  zeros(T,n+1),    # dsndb
  zeros(T,n+1,n),  # dcdu
  zeros(T,n+3),    # dfdrbc
  zeros(T,n+2),    # dfdrbu
  nmax,		   # nmax
  zeros(T,nmax),   # Iv_coeff
  zeros(T,2,2,nmax),   # Jv_coeff for k^2 < 1 & k^2 > 1; v_max & v_max-1; series coefficients
  zeros(T,2,2,nmax),   # dJvdk_coeff for k^2 < 1 & k^2 > 1; v_max & v_max-1; series coefficients
  zero(T),         # k^2
  zero(T),         # k
  zero(T),         # k_c
  zero(T),         # k_c*k
  zero(T),         # kappa/kap0
  zero(T),         # pimkap1/pi-kappa_1
  zero(T),         # kite_area2/2*A_kite
  zero(T),         # E(k^2)
  zero(T),         # (E(m)-(1-m)K(m))/m
  zero(T),         # 1-(b-r)^2
  zero(T),         # 1-(b+r)^2
  zero(T),         # 1/(1-(b-r)^2)
  zero(T),         # sqrt(1-(b-r)^2)
  zero(T),         # 4*b*r
  zero(T),	   # sqrt(b*r)
  zero(T),         # 1/sqrt(b*r)
  zero(T),         # 1/(4*b*r)
  zero(T),         # 1/k^2
  zero(T),	   # 1-k^2 (or 1-1/k^2 if k > 1)
  zero(T),         # sqrt(1-r^2)
  zero(T),         # den = 1/(pi*(c[1] + c[2]*2/3))
  one(T)/3,        # 1/3
  convert(T,2)/3,  # 2/3
  zero(T),         # sqrt(r*(1-r))
  zeros(v_max+1,v_max) # binomial(n,i)
)
# Initialize the series coefficients for I_{v_max}:
Iv_series_coeff!(trans)
# Initialize the series coefficients for J_{v_max} and J_{v_max-1}, and if
# t.grad is true, will also compute coefficients for dJ/dk_{v_max} and _{v_max-1}:
dJvdk_series_coeff!(trans)
# Compute binomial coefficients:
for n=1:trans.v_max
  trans.bincoeff[1,n]=one(T)
  trans.bincoeff[n+1,n]=one(T)
  for i=1:n-1
    if n == 1
      trans.bincoeff[i+1,n]=binomial(n,i)
    else
      trans.bincoeff[i+1,n]=trans.bincoeff[i,n-1]+trans.bincoeff[i+1,n-1]
    end
  end
end  
if grad
  compute_c_n_grad!(trans)
else
  compute_c_n!(trans)
end
trans.den = inv(pi*(trans.c_n[1]+trans.twothird*trans.c_n[2]))  # for c_2 and above, the flux is zero.
return trans
end
