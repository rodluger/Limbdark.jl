function Mm_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_m. (don't need its derivative, dM_m/dk.)
# Computing leading coefficient (n=0).  Need to compute for
# m_max-3 to m_max.  Also, need to keep track of k^2 < 1
# and k^2 >= 1.  So, Mm_coeff is 3-dimensional array:
# Mm_coeff[2,4,nmax], where first dimension is for k^2 < 1 (1)
# and k^2 >= 1 (2); second dimension is for m_max (1) to m_max-3 (4),
# while third dimension is for the series coefficients.
# Computing leading coefficient (n=0):
for k2 = 0:1
coeff = zero(T)
# Loop over m_max to m_max-3:
for j=1:4
  # m is m_max to m_max-3 (we need both since downward recursion requires two):
  m = t.m_max-(j-1)
  mhalf = 0.5*m
  if k2 < 1
    coeff = sqrt(pi)*exp(lgamma(mhalf+1.0)-lgamma(mhalf+1.5))
    # Add leading term to M_m:
    t.Mm_coeff[1,j,1] = coeff
#    if t.grad
#      t.dMmdk_coeff[1,j,1] = coeff*(m+1)
#    end
    # Now, compute higher order terms until desired precision is reached:
    @inbounds for i=1:t.nmax-1
#      coeff *= (1.5-i)^2/(mhalf+2.5-i)
      coeff *= (2*i-1)^2/(2*i*(1+m+2*i))
      t.Mm_coeff[1,j,i+1] = coeff
#      if t.grad
#        t.dMmdk_coeff[1,j,i+1] = coeff*(2j+m+1)
#      end
    end
  else # k^2 >= 1
    coeff = convert(T,pi)
    # Store leading terms:
    t.Mm_coeff[2,j,1] = coeff
#    if t.grad
#      t.dMmdk_coeff[2,j,1] = coeff*m
#    end
    # Loop over higher order terms:
    if iseven(m)
      # If m is even, then series truncates:
      jmax = div(m,2)
    else
      jmax = t.nmax-1
    end
    for i = 1:jmax
#      coeff *= (1.5-i)*(-mhalf-i+1.0)/(2.0-i)
      coeff *= (2+m-2*i)*(2*i-1)/(4*i^2)
      t.Mm_coeff[2,j,i+1] = coeff
#      if t.grad
#        t.dMmdk_coeff[2,j,i+1] = -2j*coeff
#      end
    end
  end
end
end
end
