function compute_c_n(u_n::Array{T,1}) where {T <: Real}
# Transform the u_n coefficients to c_n, which are coefficients
# of the basis in which the P(G_n) functions are computed.
n = length(u_n)
c_n = zeros(T,n+3)
a_n = zeros(T,n+1)
# Pre-allocated memory - just need to zero terms:
a_n[1] = one(T)  # Add in the first constant coefficient term
for i=1:n
  # Compute the contribution to a_n*\mu^n
  bcoeff = one(T)
  for j=0:i
#    a_n[j+1] -= u_n[i]*binomial(i,j)*(-1)^j
    a_n[j+1] -= u_n[i]*bcoeff*(-1)^j
    bcoeff *= (i-j)/(j+1)
#    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
  end
end
# Now, compute the c_n coefficients:
for j=n:-1:2
  c_n[j+1] = a_n[j+1]/(j+2)+c_n[j+3]
end
c_n[2] = a_n[2]+3*c_n[4]
c_n[1] = a_n[1]+2*c_n[3]
return c_n[1:n+1]
end

function compute_c_n_grad(u_n::Array{T,1}) where {T <: Real}
# Transform the u_n coefficients to c_n, which are coefficients
# of the basis in which the P(G_n) functions are computed.
# Compute the derivatives of the flux with respect to the u coefficients.
n = length(u_n)
# We define c_n with two extra elements which are zero:
c_n = zeros(T,n+3)
a_n = zeros(T,n+1)
dadu = zeros(T,n+1,n)
dcdu = zeros(T,n+3,n)
a_n[1] = one(T)  # Add in the first constant coefficient term
for i=1:n
  # Compute the contribution to a_n*\mu^n
  bcoeff = one(T)
  for j=0:i
#    a_n[j+1] -= u_n[i]*binomial(i,j)*(-1)^j
    a_n[j+1] -= u_n[i]*bcoeff*(-1)^j
#    dadu[j+1,i] -= binomial(i,j)*(-1)^j
    dadu[j+1,i] -= bcoeff*(-1)^j
    bcoeff *= (i-j)/(j+1)
#    println("i: ",i," j: ",j," a_i: ",a_n[j+1])
  end
end
# Now, compute the c_n coefficients and propagate derivatives:
for j=n:-1:2
  c_n[j+1] = a_n[j+1]/(j+2)+c_n[j+3]
  for i=1:n
    dcdu[j+1,i] = dadu[j+1,i]/(j+2) + dcdu[j+3,i]
  end
end
c_n[2] = a_n[2]+3*c_n[4]
for i=1:n
  dcdu[2,i] = dadu[2,i] + 3*dcdu[4,i]
end
c_n[1] = a_n[1]+2*c_n[3]
for i=1:n
  dcdu[1,i] = dadu[1,i] + 2*dcdu[3,i]
end
return c_n[1:n+1],dcdu[1:n+1,:]
end
