# Series evaluation of:
#    M_m/(4*b*r)^{m/2} = \int_{-\kappa/2}^{\kappa/2} (k^2-\sin^2{x})^{m/2} dx
# for the top four values of m = m_max-3 to m_max.
function Mm_series!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_m:
k2 = t.k2; m_max = t.m_max
# Computing leading coefficient (n=0):
if k2 < 1
  tol = eps(k2); term = zero(T)
  fac = t.sqonembmr2^(m_max-3)*t.k
# Now, compute higher order terms until desired precision is reached:
  @inbounds for j=1:4
    # Add leading term to M_m:
    t.Mm[m_max-3+j] = t.Mm_coeff[1,j,1]
    k2n = one(T)   # k^{2n}
    @inbounds for n =1:t.nmax-1
      k2n *= k2
      term = k2n*t.Mm_coeff[1,j,n+1]
      t.Mm[m_max-3+j] += term
      if abs(term) < tol
        break
      end
    end
    t.Mm[m_max-3+j] *= fac
    fac *= t.sqonembmr2
  end
  return
else # k^2 >= 1
  t.k2inv = inv(k2)
  tol = eps(t.k2inv); term = zero(T)
  fac = t.sqonembmr2^(m_max-3)
  @inbounds for j=1:4
    t.Mm[m_max-3+j] = t.Mm_coeff[2,j,1]
    k2n = one(T)
    @inbounds for n = 1:t.nmax-1
      k2n *= t.k2inv
      term = k2n*t.Mm_coeff[2,j,n+1]
      t.Mm[m_max-3+j] += term
      if abs(term) < tol
        break
      end
    end
    t.Mm[m_max-3+j] *= fac
    fac *= t.sqonembmr2
  end
  return
end
end

# Recursive computation of M_m starting at m=1 to 4, and raising to m_max:
function Mm_raise!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; m_max=t.m_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2); sqarea=t.sqarea; onemr2mb2=t.onemr2mb2
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first four integrals:
Mm_four!(t)
@inbounds for m=4:m_max
  t.Mm[m+1]=(2*(m-1)*onemr2mb2*t.Mm[m-1]+(m-2)*sqarea*t.Mm[m-3])*t.minv[m]
end
return
end

# Compute M_m for m=0 to 3:
function Mm_four!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; m_max=t.m_max; k2 = t.k2; kc = t.kc; t.k2inv = inv(t.k2)
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute first four integrals:
if k2 < 1.0
  t.Mm[1] = t.kap0
  t.Mm[2] = 2*t.sqbr*2*k2*t.Em1mKdm
  t.Mm[3] = t.Mm[1]*t.onemr2mb2+t.kite_area2
  t.Mm[4] = (2*t.sqbr)^3*t.twothird*(k2*t.Eofk+k2*(3*k2-2)*t.Em1mKdm)
else
  t.Mm[1] = pi
  t.Mm[2] = 2*t.sqonembmr2*t.Eofk
  t.Mm[3] = pi*t.onemr2mb2
  t.Mm[4] = t.sqonembmr2^3*t.twothird*((3-2*t.k2inv)*t.Eofk+t.k2inv*t.Em1mKdm)
end
return
end

# Recursive computation of M_m starting at m=m_max to m_max-3, and lowering to m=1:
function Mm_lower!(t::Transit_Struct{T}) where {T <: Real}
r=t.r; b=t.b; m_max=t.m_max; k2 = t.k2
# Computes the integrals:
# M_m(r,b) = (4*b*r)^(m/2) \int_{-\kappa/2}^{\kappa/2} (k^2 - \sin^2{x})^{m/2} dx
# where k^2 = (1-(r-b)^2)/(4*b*r), \kappa/2 = \sin^{-1}(k) for k<1; otherwise \kappa=\pi.
# Compute series version:
Mm_series!(t)
# Now iterate downwards:
invsqarea = inv(t.sqarea)
@inbounds for m=m_max-4:-1:4
  t.Mm[m+1]=((m+4)*t.Mm[m+5]-2*(m+3)*t.onemr2mb2*t.Mm[m+3])*invsqarea*t.minv[m+2]
end
# Now, compute lowest four exactly:
Mm_four!(t)
return
end
