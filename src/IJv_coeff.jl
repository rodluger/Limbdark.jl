# Computes coefficients for the I_v and J_v series expansions.

function Iv_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute coefficients in the
# series solution of I_{v_max}, needed for k^2 < 1:
v = t.v_max
# Computing leading coefficient (n=0):
coeff = 2/(2v+one(T))
# Compute leading coefficient in I_v series:
t.Iv_coeff[1] = convert(T,coeff)
# Now, compute higher order terms until desired precision is reached:
for i = 1:t.nmax-1
  n = 2i
  coeff *= (n-1)*(n+2v-1)
  coeff /= n*(n+2v+1)
  t.Iv_coeff[i+1] = coeff
end
return
end

function dJvdk_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute J_v and its derivative, dJ_v/dk.
# Computing leading coefficient (n=0).  Need to compute for
# both v_max and v_max-1.  Also, need to keep track of k^2 < 1
# and k^2 >= 1.  So, Jv_coeff is 3-dimensional array:
# Jv_coeff[2,2,nmax], where first dimension is for k^2 <1 (1)
# and k^2 >= 1 (2); second dimension is for v_max (1) and v_max-1 (2),
# while third dimension is for the series coefficients.
# Computing leading coefficient (n=0):
k2 = t.k2
coeff = zero(T)
# Loop over v_max and v_max-1:
for j=1:2
  # v is v_max or v_max-1 (we need both since downward recursion requires two):
  v = t.v_max-(j-1)
  if k2 < 1
    coeff = 0.75*convert(T,pi)/exp(lfactorial(v+2))
    # multiply by (2v-1)!!
    @inbounds for i=2:2:2v
      coeff *= 0.5*(i-1)
    end
    # Add leading term to J_v:
    t.Jv_coeff[1,j,1] = coeff
    if t.grad
      t.dJvdk_coeff[1,j,1] = coeff*(2v+1)
    end
    # Now, compute higher order terms until desired precision is reached:
    @inbounds for i=1:t.nmax-1
      n = 2i
      coeff *= (n-1)*(n+2v-1)
      coeff /= n*(n+2v+4)
      t.Jv_coeff[1,j,i+1] = coeff
      if t.grad
        t.dJvdk_coeff[1,j,i+1] = coeff*(n+2v+1)
      end
    end
  else # k^2 >= 1
    coeff = convert(T,pi)
    # Compute (2v-1)!!/(2^v v!):
    @inbounds for i=2:2:2v
      coeff *= (i-1)/i
    end
    # Store leading terms:
    t.Jv_coeff[2,j,1] = coeff
    if t.grad
      t.dJvdk_coeff[2,j,1] = zero(T)
    end
    # Loop over higher order terms:
    for i = 1:t.nmax-1
      n = 2i
      coeff *= (n-5)*(n+2v-1)
      coeff /= (n*(n+2v))
      t.Jv_coeff[2,j,i+1] = coeff
      if t.grad
        t.dJvdk_coeff[2,j,i+1] = -n*coeff
      end
    end
  end
end
end
