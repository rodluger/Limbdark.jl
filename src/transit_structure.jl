
"""
    Transit_Struct

Structure to hold arrays and other quantities for computing transit. 

# Members
- `r::T`: radius ratio
- `b::T`: impact parameter
- `u_n::Array{T,1}`: limb-darkening coefficients
- `n::Int64`: number of limb-darkening coefficients
- `n_max::Int64`: maximum number of terms in M_n
- `g_n::Array{T,1}`: Green's basis coefficients
- `sn::Array{T,1}`: Green's basis terms
- `Mn::Array{T,1}`: integral over ``(k^2-\\sin^2(\\phi))^{m/2}``
- `Nn::Array{T,1}`: integral over ``(k^2-\\sin^2(\\phi))^{m/2}\\sin^2(\\phi)``
- `grad::Bool`: `true` for gradient; `false` for no gradient
- `s2_grad::Array{T,1}`: gradient of ``s2`` with respect to ``(r,b)`` - this is ``S_1 = s_2`` from starry
- `dsndr::Array{T,1}`: derivatives of Green's basis with respect to ``r``
- `dsndb::Array{T,1}`: derivatives of Green's basis with respect to ``b``
- `dgdu::Array{T,2}`: derivatives ``g_n`` with respect to ``u_n``
- `dfdrb::Array{T,1}`: derivative of flux with respect to `r,b`
- `dfdg::Array{T,1}`: derivative of flux with respect to ``g_n``
- `dfdu::Array{T,1}`: derivative of flux with respect to `u_n`
- `jmax::Int64`: maximum number of terms in series expansions of ``I_v`` and ``J_v``
- `Mn_coeff::Array{T,3}`: coefficients for series expansion of ``M_n``
- `Nn_coeff::Array{T,2}`: coefficients for series expansion of ``N_n``
- `ninv::Array{T,1}`: inverse of the integers ``n``
- `k2::T`: ``k^2 = (1-(r-b)^2)/(4br)``
- `k::T`: ``k = \\sqrt{k^2}``
- `kc::T`: ``k_c = \\sqrt{1-k^2}`` (unless ``k > 1``, then it is ``\\sqrt{1-1/k^2}``)
- `kck::T`: ``k_c k``
- `kap0::T`: ``\\kappa = \\sin^{-1}(k)`` (``=\\kappa_0`` in M&A)
- `pimkap1::T`: ``\\pi-\\kappa_1``
- `sqarea::T`: ``(1-(b-r)^2)((b+r)^2-1)``
- `kite_area2::T`: ``2A_{kite}``
- `Eofk::T`: ``E(k^2)`` is complete elliptic integral of first kind
- `Em1mKdm::T`: ``(E(m)-(1-m)K(m))/m`` is complete elliptic integral with ``m=k^2``
- `onembmr2:: T`: ``1-(b-r)^2``
- `onembpr2:: T`: ``1-(b+r)^2``
- `onembmr2inv::T`: ``1/(1-(b-r)^2)``
- `onemr2mb2::T`: ``1-r^2-b^2``
- `sqonembmr2::T`: ``\\sqrt{1-(b-r)^2}``
- `fourbr::T`: ``4br``
- `sqbr::T`: ``\\sqrt{br}``
- `sqbrinv::T`: ``1/\\sqrt{br}``
- `fourbrinv::T`: ``1/(4br)``
- `k2inv::T`: ``1/k^2``
- `kc2::T`: ``1-k^2`` (or ``1-1/k^2`` if ``k > 1``)
- `sqrt1mr2:: T`: ``\\sqrt{1-r^2}``
- `den::T`: ``1/(g_1 + 2/3 g_2)``
- `third::T`: ``1/3``
- `twothird:: T`: ``2/3``
- `sqr1mr::T`: ``\\sqrt{r(1-r)}`` if ``r < 1``
"""
mutable struct Transit_Struct{T}
  # Structure to hold arrays and other quantities for computing transit:
  r      ::T           # radius ratio
  b      ::T           # impact parameter
  u_n    ::Array{T,1}  # limb-darkening coefficients
  n      ::Int64       # number of limb-darkening coefficients
  n_max  ::Int64       # maximum number of terms in M_n
  g_n    ::Array{T,1}  # Green's basis coefficients
  sn     ::Array{T,1}  # Green's basis terms
  Mn     :: Array{T,1} # integral over (k^2-sin^2(phi))^(m/2)
  Nn     :: Array{T,1} # integral over (k^2-sin^2(phi))^(m/2)*sin^2(phi)
  grad   ::Bool        # true for gradient; false for no gradient
  s2_grad::Array{T,1}  # gradient of s2 with respect to (r,b) - this is S_1 = s_2 from starry
  dsndr  :: Array{T,1} # derivatives of Green's basis with respect to r
  dsndb  :: Array{T,1} # derivatives of Green's basis with respect to b
  dgdu   ::Array{T,2}  # derivatives g_n with respect to u_n
  dfdrb  ::Array{T,1}  # derivative of flux with respect to r,b
  dfdg   ::Array{T,1}  # derivative of flux with respect to g_n
  dfdu   ::Array{T,1}  # derivative of flux with respect to u_n
  jmax   ::Int64       # maximum number of terms in series expansions of I_v and J_v
  Mn_coeff::Array{T,3} # coefficients for series expansion of M_n
  Nn_coeff::Array{T,2} # coefficients for series expansion of N_n
  ninv   ::Array{T,1}  # inverse of the integers n
  k2     ::T           # k^2 = (1-(r-b)^2)/(4br)
  k      ::T           # k = sqrt(k^2)
  kc     ::T           # k_c = sqrt(1-k^2) (unless k > 1, then it is sqrt(1-1/k^2))
  kck    ::T           # k_c * k
  kap0   ::T           # kappa = sin^{-1}(k) (=kappa_0 in M&A)
  pimkap1::T           # pi-kappa_1
  sqarea::T            # (1-(b-r)^2)*((b+r)^2-1)
  kite_area2::T        # 2*A_kite = 
  Eofk   ::T           # E(k^2) is complete elliptic integral of first kind
  Em1mKdm::T           # (E(m)-(1-m)K(m))/m is complete elliptic integral with m=k^2
  onembmr2:: T           # 1-(b-r)^2
  onembpr2:: T           # 1-(b+r)^2
  onembmr2inv::T       # 1/(1-(b-r)^2)
  onemr2mb2  ::T       # 1-r^2-b^2
  sqonembmr2 ::T       # sqrt(1-(b-r)^2)
  fourbr ::T           # 4*b*r
  sqbr   ::T           # sqrt(b*r)
  sqbrinv::T           # 1/sqrt(b*r)
  fourbrinv ::T           # 1/(4*b*r)
  k2inv  ::T           # 1/k^2
  kc2    ::T           # 1-k^2 (or 1-1/k^2 if k > 1)
  sqrt1mr2:: T           # sqrt(1-r^2)
  den    ::T           # 1/(g_1 + g_2*2/3)
  third  ::T           # 1/3
  twothird:: T           # 2/3
  sqr1mr ::T           # sqrt(r*(1-r)) if r < 1
end

using SpecialFunctions
include("compute_g_n_struct.jl")
include("Mn_coeff.jl")
include("Nn_coeff.jl")


"""
    transit_init(r,b,u_n,grad)

Initialize and return a transit structure. This structure holds user-facing
information about the transit as well as internal temporary variables
used for computing the light curve. Pass an instance of this structure
to `transit_poly()` to compute the actual flux.

# Arguments
- `r::Real`: The radius of the occultor in units of the radius of the occulted body.
- `b::Real`: The (initial) impact parameter of the occultation.
- `u_n::Array{Real,1}`: The array of limb darkening coefficients.
- `grad::Bool`: Compute the gradient of the flux as well?
"""
function transit_init(r::T,b::T,u_n::Array{T,1},grad::Bool) where {T <: Real}
  # Initializs a transit structure.
  n = length(u_n)
  # Maximum number of M_n integrals is n+2:
  n_max = n
  jmax = 100
  trans = Transit_Struct{T}(r,b,u_n,n,n_max,
    zeros(T,n+1),    # g_n
    zeros(T,n+1),    # sn
    zeros(T,n_max+1),# M_n from m= 0 to n_max
    zeros(T,n_max+1),# N_n from m= 0 to n_max
    grad,            # grad
    zeros(T,2),      # s2_grad
    zeros(T,n+1),    # dsndr
    zeros(T,n+1),    # dsndb
    zeros(T,n+1,n),  # dgdu
    zeros(T,2),      # dfdrb
    zeros(T,n+1),    # dfdg
    zeros(T,n),      # dfdu
    jmax,		   # jmax
    zeros(T,2,4,jmax),   # Mn_coeff for k^2 < 1 & k^2 > 1; n_max-3 to n_max; series coefficients
    zeros(T,2,jmax),  # Nn_coeff for k^2 < 1; n_max-1 to n_max; series coefficients
    zeros(T,n_max+2), # inverse of integers
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
  # Initialize series coefficients for M_n and N_n:
  Mn_series_coeff!(trans)
  Nn_series_coeff!(trans)
  # Inverse of integers:
  for n=1:n_max+2
    trans.ninv[n] = inv(n)
  end
  if grad
    compute_g_n_grad!(trans)
  else
    compute_g_n!(trans)
  end
  trans.den = inv(pi*(trans.g_n[1]+trans.twothird*trans.g_n[2]))  # for g_2 and above, the flux is zero.
  return trans
end
