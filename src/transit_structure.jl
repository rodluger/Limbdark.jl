
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
  grad    :: Bool        # true for gradient; false for no gradient
  dIvdk   ::  Array{T,1} # derivative of I_v with respect to k
  dJvdk   ::  Array{T,1} # derivative of J_v with respect to k
  s2_grad :: Array{T,1}  # gradient of s2 with respect to (r,b) - this is S_1 = s_2 from starry
  dsndr   ::  Array{T,1} # derivatives of Green's basis with respect to r
  dsndb   ::  Array{T,1} # derivatives of Green's basis with respect to b
  dcdu    :: Array{T,2}  # derivatives c_n with respect to u_n
  dfdrbc  :: Array{T,1}  # derivative of flux with respect to c_n
  dfdrbu  :: Array{T,1}  # derivative of flux with respect to u_n
  nmax    :: Int64       # maximum number of terms in series expansions of I_v and J_v
  Iv_coeff :: Array{T,1} # coefficients for series expansion of I_v
  k2      :: T           # k^2 = (1-(r-b)^2)/(4br)
  k       :: T           # k = sqrt(k^2)
  kc      :: T           # k_c = sqrt(1-k^2) (unless k > 1, then it is sqrt(1-1/k^2))
  kck     :: T           # k_c * k
  kap     :: T           # kappa = sin^{-1}(k)
  Eofk    :: T           # E(k^2) is complete elliptic integral of first kind
  Em1mKdm :: T           # (E(m)-(1-m)K(m))/m is complete elliptic integral with m=k^2
end

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
  grad,            # grad
  zeros(T,v_max+1),# dIvdk
  zeros(T,v_max+1),# dJvdk
  zeros(T,2),      # s2_grad
  zeros(T,n+1),    # dsndr
  zeros(T,n+1),    # dsndb
  zeros(T,n+1,n),  # dcdu
  zeros(T,n+3),    # dfdrbc
  zeros(T,n+2),    # dfdrbu
  nmax,		   # nmax
  zeros(T,nmax),   # Iv_coeff
  zero(T),         # k^2
  zero(T),         # k
  zero(T),         # k_c
  zero(T),         # k_c*k
  zero(T),         # kappa
  zero(T),         # E(k^2)
  zero(T)          # (E(m)-(1-m)K(m))/m
)
# Initialize the series coefficients for I_{v_max}:
Iv_series_coeff!(trans)
if grad
  compute_c_n_grad!(trans)
else
  compute_c_n!(trans)
end
return trans
end
