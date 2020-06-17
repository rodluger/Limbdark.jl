function Nn_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute N_n.
# Computing leading coefficient (n=0).  Need to compute for
# n_max-3 to n_max.  Also, need just for k^2 < 1
# Nn_coeff is 2-dimensional array:
# Nn_coeff[2,jmax], where first dimension is for n_max-1 (1) to n_max (2),
# while third dimension is for the series coefficients.
# Computing leading coefficient (n=0):
coeff = zero(T)
# Loop over n_max-1 to n_max:
for j=1:2
  # m is n_max-1 to n_max (we need both since downward recursion requires two):
  m = t.n_max+j-2
  mhalf = 0.5*m
  coeff = 0.5*sqrt(pi)*exp(logabsgamma(mhalf+1.0)[1]-logabsgamma(mhalf+2.5)[1])
  # Add leading term to N_n:
  t.Nn_coeff[j,1] = coeff
  # Now, compute higher order terms until desired precision is reached:
  @inbounds for i=1:t.jmax-1
    coeff *= (4*i^2-1)/(2*i*(3+m+2*i))
    t.Nn_coeff[j,i+1] = coeff
  end
end
end
