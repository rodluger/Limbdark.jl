mutable struct Transit_Struct{T}
  # Structure to hold arrays and other quantities for computing transit:
  r       :: T           # radius ratio
  b       :: T           # impact parameter
  u_n     :: Array{T,1}  # limb-darkening coefficients
  n       :: Int64       # number of limb-darkening coefficients
  m_max   :: Int64       # maximum number of terms in M_m
  d_n     :: Array{T,1}  # Green's basis coefficients
  sn      :: Array{T,1}  # Green's basis terms
  Mn      ::  Array{T,1} # integral over (k^2-sin(phi)^2)^(m/2)
  grad    :: Bool        # true for gradient; false for no gradient
  s2_grad :: Array{T,1}  # gradient of s2 with respect to (r,b) - this is S_1 = s_2 from starry
  dsndr   ::  Array{T,1} # derivatives of Green's basis with respect to r
  dsndb   ::  Array{T,1} # derivatives of Green's basis with respect to b
  dcdu    :: Array{T,2}  # derivatives d_n with respect to u_n
  dfdrb   :: Array{T,1}  # derivative of flux with respect to r,b
  dfdc    :: Array{T,1}  # derivative of flux with respect to d_n
  dfdu    :: Array{T,1}  # derivative of flux with respect to u_n
  nmax    :: Int64       # maximum number of terms in series expansions of I_v and J_v
  Mn_coeff :: Array{T,3} # coefficients for series expansion of M_m
  minv    :: Array{T,1}  # inverse of the integers m
  k2      :: T           # k^2 = (1-(r-b)^2)/(4br)
  k       :: T           # k = sqrt(k^2)
  kc      :: T           # k_c = sqrt(1-k^2) (unless k > 1, then it is sqrt(1-1/k^2))
  kck     :: T           # k_c * k
  kap0    :: T           # kappa = sin^{-1}(k) (=kappa_0 in M&A)
  pimkap1 :: T           # pi-kappa_1
  sqarea :: T            # (1-(b-r)^2)*((b+r)^2-1)
  kite_area2 :: T        # 2*A_kite = 
  Eofk    :: T           # E(k^2) is complete elliptic integral of first kind
  Em1mKdm :: T           # (E(m)-(1-m)K(m))/m is complete elliptic integral with m=k^2
  onembmr2:: T           # 1-(b-r)^2
  onembpr2:: T           # 1-(b+r)^2
  onembmr2inv :: T       # 1/(1-(b-r)^2)
  onemr2mb2   :: T       # 1-r^2-b^2
  sqonembmr2  :: T       # sqrt(1-(b-r)^2)
  fourbr  :: T           # 4*b*r
  sqbr    :: T           # sqrt(b*r)
  sqbrinv :: T           # 1/sqrt(b*r)
  fourbrinv  :: T           # 1/(4*b*r)
  k2inv   :: T           # 1/k^2
  kc2     :: T           # 1-k^2 (or 1-1/k^2 if k > 1)
  sqrt1mr2:: T           # sqrt(1-r^2)
  den     :: T           # 1/(d_1 + d_2*2/3)
  third   :: T           # 1/3
  twothird:: T           # 2/3
  sqr1mr  :: T           # sqrt(r*(1-r)) if r < 1
end

using SpecialFunctions
include("compute_d_n_struct.jl")
include("Mn_coeff.jl")

function transit_init(r::T,b::T,u_n::Array{T,1},grad::Bool) where {T <: Real}
# Initializs a transit structure.
n = length(u_n)
# Maximum number of M_m integrals is n+2:
m_max = n+1
nmax = 100
trans = Transit_Struct{T}(r,b,u_n,n,m_max,
  zeros(T,n+1),    # d_n
  zeros(T,n+1),    # sn
  zeros(T,m_max+1),# M_m from m= 0 to m_max
  grad,            # grad
  zeros(T,2),      # s2_grad
  zeros(T,n+1),    # dsndr
  zeros(T,n+1),    # dsndb
  zeros(T,n+1,n),  # dcdu
  zeros(T,2),      # dfdrb
  zeros(T,n+1),    # dfdc
  zeros(T,n),      # dfdu
  nmax,		   # nmax
  zeros(T,2,4,nmax),   # Mn_coeff for k^2 < 1 & k^2 > 1; m_max-3 to m_max; series coefficients
  zeros(T,m_max+1), # inverse of integers
  zero(T),         # k^2
  zero(T),         # k
  zero(T),         # k_c
  zero(T),         # k_c*k
  zero(T),         # kappa/kap0
  zero(T),         # pimkap1/pi-kappa_1
  zero(T),         # sqarea
  zero(T),         # kite_area2/2*A_kite
  zero(T),         # E(k^2)
  zero(T),         # (E(m)-(1-m)K(m))/m
  zero(T),         # 1-(b-r)^2
  zero(T),         # 1-(b+r)^2
  zero(T),         # 1/(1-(b-r)^2)
  zero(T),         # 1-r^2-b^2
  zero(T),         # sqrt(1-(b-r)^2)
  zero(T),         # 4*b*r
  zero(T),	   # sqrt(b*r)
  zero(T),         # 1/sqrt(b*r)
  zero(T),         # 1/(4*b*r)
  zero(T),         # 1/k^2
  zero(T),	   # k_c^2 = 1-k^2 (or k_c^2=1-1/k^2 if k > 1)
  zero(T),         # sqrt(1-r^2)
  zero(T),         # den = 1/(pi*(c[1] + c[2]*2/3))
  one(T)/3,        # 1/3
  convert(T,2)/3,  # 2/3
  zero(T)         # sqrt(r*(1-r))
)
# Initialize series coefficients for M_m:
Mn_series_coeff!(trans)
# Inverse of integers:
for m=1:m_max+1
  trans.minv[m] = inv(m)
end
if grad
  compute_d_n_grad!(trans)
else
  compute_d_n!(trans)
end
trans.den = inv(pi*(trans.d_n[1]+trans.twothird*trans.d_n[2]))  # for d_2 and above, the flux is zero.
return trans
end
