function Mn_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute M_n. (don't need its derivative, dM_n/dk.)
# Computing leading coefficient (n=0).  Need to compute for
# n_max-3 to n_max.  Also, need to keep track of k^2 < 1
# and k^2 >= 1.  So, Mn_coeff is 3-dimensional array:
# Mn_coeff[2,4,jmax], where first dimension is for k^2 < 1 (1)
# and k^2 >= 1 (2); second dimension is for n_max-3 (1) to n_max (4),
# while third dimension is for the series coefficients.
# Computing leading coefficient (n=0):
for k2 = 0:1
coeff = zero(T)
# Loop over n_max to n_max-3:
for j=1:4
  # m is n_max to n_max-3 (we need both since downward recursion requires two):
  m = t.n_max+j-4
  mhalf = 0.5*m
  if k2 < 1
      coeff = sqrt(pi)*exp(logabsgamma(mhalf+1.0)[1]-logabsgamma(mhalf+1.5)[1])
    # Add leading term to M_n:
    t.Mn_coeff[1,j,1] = coeff
    # Now, compute higher order terms until desired precision is reached:
    @inbounds for i=1:t.jmax-1
      coeff *= (2*i-1)^2/(2*i*(1+m+2*i))
      t.Mn_coeff[1,j,i+1] = coeff
    end
  else # k^2 >= 1
    coeff = convert(T,pi)
    # Store leading terms:
    t.Mn_coeff[2,j,1] = coeff
    # Loop over higher order terms:
    if iseven(m)
      # If m is even, then series truncates:
      jmax = div(m,2)
    else
      jmax = t.jmax-1
    end
    for i = 1:jmax
      coeff *= (2+m-2*i)*(1-2*i)/(4*i^2)
      t.Mn_coeff[2,j,i+1] = coeff
    end
  end
end
end
end
