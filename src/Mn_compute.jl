
# Series evaluation of:
#    M_n/(4*b*r)^{n/2} = \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2{x})^{n/2} dx
# for the top four values of n = n_max-3 to n_max.
function Mn_series!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_n:
k2 = t.k2; n_max = t.n_max
# Computing leading coefficient (n=0):
if k2 < 1
  tol = eps(k2); term = zero(T)
  fac = t.sqonembmr2^(n_max-3)*t.k
# Now, compute higher order terms until desired precision is reached:
  @inbounds for j=1:4
    # Add leading term to M_n:
    t.Mn[n_max-3+j] = t.Mn_coeff[1,j,1]
    k2n = one(T)   # k^{2n}
    @inbounds for n =1:t.jmax-1
      k2n *= k2
      term = k2n*t.Mn_coeff[1,j,n+1]
      t.Mn[n_max-3+j] += term
      if abs(term) < tol
        break
      end
    end
    t.Mn[n_max-3+j] *= fac
    fac *= t.sqonembmr2
  end
  return
else # k^2 >= 1
  t.k2inv = inv(k2)
  tol = eps(t.k2inv); term = zero(T)
  fac = t.sqonembmr2^(n_max-3)
  @inbounds for j=1:4
    t.Mn[n_max-3+j] = t.Mn_coeff[2,j,1]
    k2n = one(T)
    @inbounds for n = 1:t.jmax-1
      k2n *= t.k2inv
      term = k2n*t.Mn_coeff[2,j,n+1]
      t.Mn[n_max-3+j] += term
      if abs(term) < tol
        break
      end
    end
    t.Mn[n_max-3+j] *= fac
    fac *= t.sqonembmr2
  end
  return
end
end

# Recursive computation of M_n starting at m=1 to 4, and raising to n_max:
function Mn_raise!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2); sqarea=t.sqarea; onemr2mb2=t.onemr2mb2
# Computes the integrals:
# M_n(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first four integrals:
Mn_four!(t)
@inbounds for m=4:n_max
  t.Mn[m+1]=(2*(m-1)*onemr2mb2*t.Mn[m-1]+(m-2)*sqarea*t.Mn[m-3])*t.ninv[m]
end
return
end

# Compute M_n for m=0 to 3:
function Mn_four!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2)
# Computes the integrals:
# M_n(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first four integrals:
if k2 < 1.0
  t.Mn[1] = t.kap0
  t.Mn[2] = 2*t.sqbr*2*k2*t.Em1mKdm
  t.Mn[3] = t.Mn[1]*t.onemr2mb2+t.kite_area2
  t.Mn[4] = (2*t.sqbr)^3*t.twothird*k2*(t.Eofk+(3*k2-2)*t.Em1mKdm)
else
  t.Mn[1] = pi
  t.Mn[2] = 2*t.sqonembmr2*t.Eofk
  t.Mn[3] = pi*t.onemr2mb2
  t.Mn[4] = t.sqonembmr2^3*t.twothird*((3-2*t.k2inv)*t.Eofk+t.k2inv*t.Em1mKdm)
end
return
end

# Recursive computation of M_n starting at m=n_max to n_max-3, and lowering to m=1:
function Mn_lower!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; n_max=t.n_max; k2 = t.k2
# Computes the integrals:
# M_n(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute series version:
Mn_series!(t)
# Now iterate downwards:
invsqarea = inv(t.sqarea)
@inbounds for m=n_max-4:-1:4
  t.Mn[m+1]=((m+4)*t.Mn[m+5]-2*(m+3)*t.onemr2mb2*t.Mn[m+3])*invsqarea*t.ninv[m+2]
end
# Now, compute lowest four exactly:
Mn_four!(t)
return
end
