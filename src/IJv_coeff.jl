# Computes coefficients for the I_v and J_v series expansions.

function Iv_series_coeff!(t::Transit_Struct{T}) where {T <: Real}
# Use series expansion to compute coefficients in the
# series solution of I_{v_max}, needed for k^2 < 1:
v = t.v_max
# Computing leading coefficient (n=0):
coeff = 2.0/(2v+1.0)
# Compute leading coefficient in I_v series:
t.Iv_coeff[1] = coeff
# Now, compute higher order terms until desired precision is reached:
for i = 1:t.nmax-1
  n = 2i
  coeff *= (n-1)*(n+2v-1)
  coeff /= n*(n+2v+1)
  t.Iv_coeff[i+1] = coeff
end
return
end

function Jv_series_coeff(k2::T,v::Int64) where {T <: Real}
# Use series expansion to compute J_v:
nmax = 100
n = 2; error = Inf; if k2 < 1; tol = eps(k2); else; tol = eps(inv(k2)); end
# Computing leading coefficient (n=0):
#coeff = 3pi/(2^(2+v)*factorial(v+2))
if k2 < 1
#  coeff = 0.75*pi/exp(lfact(v+2))
  coeff = 0.75*pi/exp(lfactorial(v+2))
# multiply by (2v-1)!!
  @inbounds for i=2:2:2v
    coeff *= (i-1)/2
  end
# Add leading term to J_v:
  Jv = convert(T,coeff)
# Now, compute higher order terms until desired precision is reached:
#  while n < nmax && abs(error) > tol
  for n=2:2:nmax
    coeff *= (n-1)*(n+2v-1)
    coeff /= n*(n+2v+4)
    coeff *= k2
    Jv += coeff
#    error = coeff
    if abs(coeff) < tol
      break
    end
#    n += 2
  end
  return Jv*k2^v*sqrt(k2)
else # k^2 >= 1
  coeff = convert(typeof(k2),pi)
  # Compute (2v-1)!!/(2^v v!):
  @inbounds for i=2:2:2v
    coeff *= (i-1)/i
  end
  Jv = convert(typeof(k2),coeff)
  k2inv = inv(k2)
#  while n < nmax && abs(error) > tol
  @inbounds for n =2:2:nmax
#    coeff *= (1.-2.5/n)*(1.-.5/(n+v))/k2
#    coeff *= (1-5/(n))*(1-1/(n+2v))/k2
#    coeff *= (1.0-5.0/float(n))*(1.0-1.0/float(n+2v))/k2
    coeff *= (n-5)*(n+2v-1)
    coeff /= k2*(n*(n+2v))  # This line takes about 27% of run time!
#    coeff /= (n*(n+2v))  # This line takes about 27% of run time!
#    coeff *= k2inv
    Jv += coeff
    if abs(coeff) < tol
      break
    end
#    error = coeff
#    n += 2
  end
  return Jv
end
end

function dJv_seriesdk_coeff(k2::T,v::Int64) where {T <: Real}
# Use series expansion to compute J_v:
nmax = 100
n = 2; error = Inf; if k2 < 1; tol = eps(k2); else; tol = eps(inv(k2)); end
# Computing leading coefficient (n=0):
#coeff = 3pi/(2^(2+v)*factorial(v+2))
coeff = zero(k2)
if k2 < 1
#  coeff = 3pi/(2^(2+v)*exp(lfact(v+2)))
  coeff = 3pi/(2^(2+v)*exp(lfactorial(v+2)))
#  println("coefficient: ",coeff)
# multiply by (2v-1)!!
  @inbounds for i=2:2:2v
    coeff *= i-1
  end
# Add leading term to J_v:
  Jv = one(k2)*coeff
  dJvdk = one(k2)*coeff*(2v+1)
# Now, compute higher order terms until desired precision is reached:
#  while n < nmax && abs(error) > tol
  @inbounds for n=2:2:nmax
    coeff *= (n-1)*(n+2v-1)
    coeff /= n*(n+2v+4)
    coeff *= k2
    Jv += coeff
    dJvdk += coeff*(n+2v+1)
#    error = coeff/Jv
#    error = coeff
    if abs(coeff) < tol
      break
    end
#    n += 2
  end
  dJvdk *= k2^v
  Jv *= k2^v*sqrt(k2)
#  println("Jv: ",Jv," dJv/dk: ",dJvdk)
  return Jv,dJvdk
else # k^2 >= 1
  coeff = convert(typeof(k2),pi)
  # Compute (2v-1)!!/(2^v v!):
  @inbounds for i=2:2:2v
    coeff *= (i-1)/i
  end
  Jv = one(k2)*coeff
  dJvdk = zero(k2)
  k2inv = inv(k2)
#  while n < nmax && abs(error) > tol
  for n = 2:2:nmax
#    coeff *= (1-5/n)*(1-1/(n+2v))*k2inv
    coeff *= (n-5)*(n+2v-1)
    coeff /= (n*(n+2v))
    coeff *= k2inv

    Jv += coeff
    dJvdk -= n*coeff
#    error = coeff/Jv
#    error = coeff
    if abs(coeff) < tol
      break
    end
#    n += 2
  end
  dJvdk /= sqrt(k2)
  return Jv,dJvdk
end
end
