function compute_d_n!(t::Transit_Struct{T}) where {T <: Real}
# Transform the u_n coefficients to d_n, which are coefficients
# of the new basis in which the P(G_n) functions are computed:
#
# \sum_{i=0}^n d_n {z^n - 3/(n+2)z}
# 
t.d_n[1] = one(T)  # Add in the first constant coefficient term
for i=1:t.n
  # Compute the contribution to d_n*\mu^n
  bcoeff = one(T)
  for j=0:i
    t.d_n[j+1] -= t.u_n[i]*bcoeff*(-1)^j
    bcoeff *= (i-j)/(j+1)
  end
end
# Now, add in odd terms:
if t.n >= 3
  for j=3:2:t.n
    t.d_n[2] += 3/(2+j)*t.d_n[j+1]
  end
end
return
end

function compute_c_n_grad!(t::Transit_Struct{T}) where {T <: Real}
# Transform the u_n coefficients to c_n, which are coefficients
# of the basis in which the P(G_n) functions are computed.
# Compute the derivatives of the flux with respect to the u coefficients.
# We define c_n with two extra elements which are zero:
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
# Now, compute the c_n coefficients and propagate derivatives:
for j=t.n:-1:2
  if j >= t.n-1
    t.c_n[j+1] = a_n[j+1]/(j+2)
    for i=1:t.n
      t.dcdu[j+1,i] = dadu[j+1,i]/(j+2)
    end
  else
    t.c_n[j+1] = a_n[j+1]/(j+2)+t.c_n[j+3]
    for i=1:t.n
      t.dcdu[j+1,i] = dadu[j+1,i]/(j+2) + t.dcdu[j+3,i]
    end
  end
end
if t.n >= 3
  t.c_n[2] = a_n[2]+3*t.c_n[4]
  for i=1:t.n
    t.dcdu[2,i] = dadu[2,i] + 3*t.dcdu[4,i]
  end
else
  t.c_n[2] = a_n[2]
  for i=1:t.n
    t.dcdu[2,i] = dadu[2,i]
  end
end
if t.n >= 2
  t.c_n[1] = a_n[1]+2*t.c_n[3]
  for i=1:t.n
    t.dcdu[1,i] = dadu[1,i] + 2*t.dcdu[3,i]
  end
else
  t.c_n[1] = a_n[1]
  for i=1:t.n
    t.dcdu[1,i] = dadu[1,i]
  end
end
return
end
