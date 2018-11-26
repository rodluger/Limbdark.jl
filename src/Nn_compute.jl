# Series evaluation of:
#    N_n/(4*b*r)^{n/2} = \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2{x})^{n/2} sin^2{x} dx
# for the top two values of n = n_max-1 to n_max.
function Nn_series!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_n:
k2 = t.k2; n_max = t.n_max
# Computing leading coefficient (n=0):
if k2 < 1
  tol = eps(k2); term = zero(T)
  fac = t.sqonembmr2^(n_max-1)*t.k*t.k2
# Now, compute higher order terms until desired precision is reached:
  @inbounds for j=1:2
    # Add leading term to M_n:
    t.Nn[n_max-1+j] = t.Nn_coeff[j,1]
    k2n = one(T)   # k^{2n}
    @inbounds for n =1:t.jmax-1
      k2n *= k2
      term = k2n*t.Nn_coeff[j,n+1]
      t.Nn[n_max-1+j] += term
      if abs(term) < tol
        break
      end
    end
    t.Nn[n_max-1+j] *= fac
    fac *= t.sqonembmr2
  end
  return
else # k^2 >= 1
  # No need for this one.
end
end

# Recursive computation of N_n starting at m=0 to 1, and raising to n_max:
function Nn_raise!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2); sqarea=t.sqarea; onemr2mb2=t.onemr2mb2
# Computes the integrals:
# N_n(r,b) = (4*b*r)^(n/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{n/2} \sin^2{x} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first two integrals:
Nn_two!(t)
@inbounds for m=2:n_max
  t.Nn[m+1]=(t.Mn[m+1]+m*t.onembpr2*t.Nn[m-1])*t.ninv[m+2]
end
return
end

# Compute N_n for m=0 to 1:
function Nn_two!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2)
# Computes the integrals:
# N_n(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} \sin^2{x} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first two integrals:
if k2 <= 1.0
  t.Nn[1] = 0.5*t.kap0-t.k*t.kc 
  t.Nn[2] = t.twothird*2*t.sqbr*t.k2*(-t.Eofk+2*t.Em1mKdm)
else
  t.Nn[1] = 0.5*pi
  t.Nn[2] =  t.twothird*2*t.sqbr*t.k*(2*t.Eofk - t.Em1mKdm)
end
return
end

# Recursive computation of N_n starting at m=n_max to n_max-1, and lowering to m=1:
function Nn_lower!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2
# Computes the integrals:
# N_n(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} \sin^2{x} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute series version for top two:
Nn_series!(t)
# Now iterate downwards:
invonembpr2 = inv(t.onembpr2)
@inbounds for m=n_max-2:-1:2
  t.Nn[m+1]=((m+4)*t.Nn[m+3]-t.Mn[m+3])*invonembpr2*t.ninv[m+2]
end
# Now, compute lowest two exactly:
Nn_two!(t)
return
end
