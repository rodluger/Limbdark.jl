function compute_g_n!(t::Transit_Struct{T}) where {T <: Real}
# Transform the u_n coefficients to g_n, which are coefficients
# of the basis in which the P(G_n) functions are computed.
a_n = zeros(T,t.n+1)
# Pre-allocated memory - just need to zero terms:
a_n[1] = one(T)  # Add in the first constant coefficient term
for i=1:t.n
  # Compute the contribution to a_n*\mu^n
  bcoeff = one(T)
  for j=0:i
#    a_n[j+1] -= t.u_n[i]*binomial(i,j)*(-1)^j
    a_n[j+1] -= t.u_n[i]*bcoeff*(-1)^j
    bcoeff *= (i-j)/(j+1)
#    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
  end
end
# Now, compute the g_n coefficients:
for j=t.n:-1:2
  if j >= t.n-1
    t.g_n[j+1] = a_n[j+1]/(j+2)
  else
    t.g_n[j+1] = a_n[j+1]/(j+2)+t.g_n[j+3]
  end
end
if t.n >= 3
  t.g_n[2] = a_n[2]+3*t.g_n[4]
else
  t.g_n[2] = a_n[2]
end
if t.n >= 2
  t.g_n[1] = a_n[1]+2*t.g_n[3]
else
  t.g_n[1] = a_n[1]
end
return
end

function compute_g_n_grad!(t::Transit_Struct{T}) where {T <: Real}
# Transform the u_n coefficients to g_n, which are coefficients
# of the basis in which the P(G_n) functions are computed.
# Compute the derivatives of the flux with respect to the u coefficients.
# We define g_n with two extra elements which are zero:
a_n = zeros(T,t.n+1)
dadu = zeros(T,t.n+1,t.n)
a_n[1] = one(T)  # Add in the first constant coefficient term
for i=1:t.n
  # Compute the contribution to a_n*\mu^n
  bcoeff = one(T)
  for j=0:i
#    a_n[j+1] -= t.u_n[i]*binomial(i,j)*(-1)^j
    a_n[j+1] -= t.u_n[i]*bcoeff*(-1)^j
#    dadu[j+1,i] -= binomial(i,j)*(-1)^j
    dadu[j+1,i] -= bcoeff*(-1)^j
    bcoeff *= (i-j)/(j+1)
#    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
  end
end
# Now, compute the g_n coefficients and propagate derivatives:
for j=t.n:-1:2
  if j >= t.n-1
    t.g_n[j+1] = a_n[j+1]/(j+2)
    for i=1:t.n
      t.dgdu[j+1,i] = dadu[j+1,i]/(j+2)
    end
  else
    t.g_n[j+1] = a_n[j+1]/(j+2)+t.g_n[j+3]
    for i=1:t.n
      t.dgdu[j+1,i] = dadu[j+1,i]/(j+2) + t.dgdu[j+3,i]
    end
  end
end
if t.n >= 3
  t.g_n[2] = a_n[2]+3*t.g_n[4]
  for i=1:t.n
    t.dgdu[2,i] = dadu[2,i] + 3*t.dgdu[4,i]
  end
else
  t.g_n[2] = a_n[2]
  for i=1:t.n
    t.dgdu[2,i] = dadu[2,i]
  end
end
if t.n >= 2
  t.g_n[1] = a_n[1]+2*t.g_n[3]
  for i=1:t.n
    t.dgdu[1,i] = dadu[1,i] + 2*t.dgdu[3,i]
  end
else
  t.g_n[1] = a_n[1]
  for i=1:t.n
    t.dgdu[1,i] = dadu[1,i]
  end
end
return
end
