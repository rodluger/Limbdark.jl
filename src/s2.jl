# Uses new formulation from limbdark paper.
include("cel_bulirsch.jl")

"""
    s2(r,b)

Compute the integral of the linear limb-darkening case.

# Examples
```julia-repl
julia> s2(0.1,0.5)
2.067294367278038
```
"""
function s2(r::T,b::T) where {T <: Real}
s_2,Eofk,Em1mKdm = s2_ell(r,b)
return s_2
end

"""
    s2_ell(r,b)
 
Compute the integral of the linear limb-darkening case, s_2 from
starry (which is S_1 in Agol & Luger) and also return two elliptic 
integrals, E(k) and (E(m)-(1-m)K(m))/m, where m = k^2.

# Examples
```julia-repl
julia> s2_ell(0.1,0.5)
(2.067294367278038, 1.4726447391114554, 0.8111640472029077)
```
"""
function s2_ell(r::T,b::T) where {T <: Real}
if T == Float64
  third = third_float
elseif T == BigFloat
  third = third_big
else
  third = convert(T,1)/convert(T,3)
end
# For now, just compute linear component:
Eofk = zero(T)
Em1mKdm = zero(T)
Lambda1 = zero(T)
if b >= 1.0+r ||  r == 0.0
  # No occultation:
  Lambda1 = zero(T)  # Case 1
elseif b <= r-1.0
  # Full occultation:
  Lambda1 = zero(T)  # Case 11
else 
  if b == 0 
    Lambda1 = -2pi*sqrt(1.0-r^2)^3 # Case 10
    Eofk = .5*pi
    Em1mKdm = .25*pi
  elseif b==r
    if r == 0.5
#      Lambda1 = pi-4/3 - 2*(b-0.5) +6*(r-0.5) # Case 6; I've added in analytic first derivaties.
      Lambda1 = pi-4*third - 2*(b-0.5) +6*(r-0.5) # Case 6; I've added in analytic first derivaties.
      Eofk = one(T)
      Em1mKdm = one(T)
    elseif r < 0.5
      m = 4r^2
      Eofk = cel_bulirsch(m,one(T),one(T),1-m)
      Em1mKdm = cel_bulirsch(m,one(T),one(T),zero(T))
#      Lambda1 = pi+2/3*((2m-3)*Eofk - m*Em1mKdm)+ # Case 5
      Lambda1 = pi+2*third*((2m-3)*Eofk - m*Em1mKdm)+ # Case 5
        (b-r)*4*r*(Eofk-2*Em1mKdm)  # Adding in first derivative
    else
      m = 4r^2; minv = inv(m)
      Eofk = cel_bulirsch(minv,one(T),one(T),1-minv)
      Em1mKdm = cel_bulirsch(minv,one(T),one(T),zero(T))
#      Lambda1 = pi+1/(3*r)*(-m*Eofk + (2m-3)*Em1mKdm) - # Case 7
      Lambda1 = pi+third/r*(-m*Eofk + (2m-3)*Em1mKdm) - # Case 7
        (b-r)*2*(2*Eofk-Em1mKdm) # Adding in first derivative
    end
  else
    onembpr2 = (1-b-r)*(1+b+r); onembmr2=(r-b+1)*(1-r+b); fourbr = 4b*r; fourbrinv = inv(fourbr)
    k2 = onembpr2*fourbrinv+1
    if (b+r) > 1.0 # k^2 < 1, Case 2, Case 8
      k2c = -onembpr2*fourbrinv; kc = sqrt(k2c)
      Piofk,Eofk,Em1mKdm = cel_bulirsch(k2,kc,(b-r)^2*k2c,zero(T),one(T),one(T),3*k2c*(b-r)*(b+r),k2c,zero(T))
#      Lambda1 = onembmr2*(Piofk+ (-3+6r^2+2*b*r)*Em1mKdm-fourbr*Eofk)/(3*sqrt(b*r))
      Lambda1 = onembmr2*(Piofk+ (-3+6r^2+2*b*r)*Em1mKdm-fourbr*Eofk)*third/sqrt(b*r)
    elseif (b+r) < 1.0  # k^2 > 1, Case 3, Case 9
      onembmr2inv=inv(onembmr2); k2inv = inv(k2); k2c =onembpr2*onembmr2inv; kc = sqrt(k2c)
      bmrdbpr = (b-r)/(b+r)
      mu = 3bmrdbpr*onembmr2inv
      p = bmrdbpr^2*onembpr2*onembmr2inv
      Piofk,Eofk,Em1mKdm = cel_bulirsch(k2inv,kc,p,1+mu,one(T),one(T),p+mu,k2c,zero(T))
#      Lambda1 = 2*sqrt(onembmr2)*(onembpr2*Piofk -(4-7r^2-b^2)*Eofk)/3
      Lambda1 = 2*sqrt(onembmr2)*(onembpr2*Piofk -(4-7r^2-b^2)*Eofk)*third
    else
      # b+r = 1 or k^2=1, Case 4 (extending r up to 1)
#      Lambda1 = 2*acos(1.0-2.0*r)-2pi*convert(T,r>.5) -(4/3*(3+2r-8r^2)+
      Lambda1 = 2*acos(1.0-2.0*r)-2pi*convert(T,r>.5) -(4*third*(3+2r-8r^2)+
          8*(r+b-1)*r)*sqrt(r*(1-r)) # Adding in first derivatives
      Eofk = one(T)
      Em1mKdm = one(T)
    end
  end
end
#flux = ((1.0-convert(T,r>b))*2pi-Lambda1)/3
flux = ((1.0-convert(T,r>b))*2pi-Lambda1)*third
return flux,Eofk,Em1mKdm
end

"""
    s2!(r,b,s2_grad)
 
Computes the linear limb-darkening case, as well as the gradient,
s2_grad=[ds_2/dr,ds_2/db] is a pre-allocated two-element array.
Returns s2 and complete elliptic integrals, E(k) and 
(E(m)-(1-m)K(m))/m, where m = k^2.

# Example
```julia-repl
julia> s2_grad=[0.0,0.0]
julia> s2!(0.1,0.5,s2_grad)
(2.067294367278038, 1.4726447391114554, 0.8111640472029077)
julia> s2_grad
2-element Array{Float64,1}:
 -0.53988  
  0.0182916
```
"""
function s2!(r::T,b::T,s2_grad::Array{T,1}) where {T <: Real}
if T == Float64
  third = third_float
elseif T == BigFloat
  third = third_big
else
  third = convert(T,1)/convert(T,3)
end
Eofk = zero(T)
Em1mKdm = zero(T)
Lambda1 = zero(T)
fill!(s2_grad,zero(T))
if b >= 1.0+r ||  r == 0.0
  # No occultation:
  Lambda1 = zero(T)  # Case 1
elseif b <= r-1.0
  # Full occultation:
  Lambda1 = zero(T)  # Case 11
else 
  if b == 0 
    sqrt1mr2 = sqrt(1.0-r^2)
    Lambda1 = -2pi*sqrt1mr2^3 # Case 10
    s2_grad[1] = -2pi*r*sqrt1mr2 # dLambda/dr (dLambda/db= 0)
    Eofk = .5*pi
    Em1mKdm = .25*pi
  elseif b==r
    if r == 0.5
#      Lambda1 = pi-4/3 # Case 6; added in analytic first derivaties.
      Lambda1 = pi-4*third # Case 6; added in analytic first derivaties.
      s2_grad[1] = -2      # dLambda/dr
#      s2_grad[2] =  2/3  # dLambda/db
      s2_grad[2] =  2*third  # dLambda/db
      Eofk = one(T)
      Em1mKdm = one(T)
    elseif r < 0.5
      m = 4r^2
      Eofk = cel_bulirsch(m,one(T),one(T),1-m)
      Em1mKdm = cel_bulirsch(m,one(T),one(T),zero(T))
#      Lambda1 = pi+2/3*((2m-3)*Eofk - m*Em1mKdm) # Case 5
      Lambda1 = pi+2*third*((2m-3)*Eofk - m*Em1mKdm) # Case 5
      s2_grad[1] = -4*r*Eofk      # Adding in first derivative dLambda/dr
#      s2_grad[2] = -4*r/3*(Eofk-2*Em1mKdm) # Adding in first derivative dLambda/db
      s2_grad[2] = -4*r*third*(Eofk-2*Em1mKdm) # Adding in first derivative dLambda/db
    else
      m = 4r^2; minv = inv(m); kc = sqrt(1.0-minv)
      Eofk = cel_bulirsch(minv,one(T),one(T),1-minv)
      Em1mKdm = cel_bulirsch(minv,one(T),one(T),zero(T))
#      Lambda1 = pi+1/(3*r)*(-m*Eofk + (2m-3)*Em1mKdm)  # Case 7
      Lambda1 = pi+third/r*(-m*Eofk + (2m-3)*Em1mKdm)  # Case 7
      s2_grad[1] = -2*Em1mKdm # dLambda/dr
#      s2_grad[2] =  2/3*(2*Eofk-Em1mKdm)  # dLambda/db
      s2_grad[2] =  2*third*(2*Eofk-Em1mKdm)  # dLambda/db
    end
  else
    onembpr2 = (1-r-b)*(1+r+b); onembmr2=(r+1-b)*(1-r+b); fourbr = 4b*r; fourbrinv=inv(fourbr)
    k2 = onembpr2*fourbrinv+1
    if (b+r) > 1.0 # k^2 < 1, Case 2, Case 8
      k2c = -onembpr2*fourbrinv; kc = sqrt(k2c); sqbr=sqrt(b*r); sqbrinv = inv(sqbr)
      Piofk,Eofk,Em1mKdm = cel_bulirsch(k2,kc,(b-r)^2*k2c,zero(T),one(T),one(T),3*k2c*(b-r)*(b+r),k2c,zero(T))
#      Lambda1 = onembmr2*(Piofk+ (-3+6r^2+2*b*r)*Em1mKdm-fourbr*Eofk)*sqbrinv/3
      Lambda1 = onembmr2*(Piofk+ (-3+6r^2+2*b*r)*Em1mKdm-fourbr*Eofk)*sqbrinv*third
      s2_grad[1] = -2r*onembmr2*Em1mKdm*sqbrinv
#      s2_grad[2] = 2r*onembmr2*(-Em1mKdm+2*Eofk)*sqbrinv/3
      s2_grad[2] = 2r*onembmr2*(-Em1mKdm+2*Eofk)*sqbrinv*third
    elseif (b+r) < 1.0  # k^2 > 1, Case 3, Case 9
      onembmr2inv = inv(onembmr2); k2inv = inv(k2); k2c =onembpr2*onembmr2inv; kc = sqrt(k2c)
      bmrdbpr = (b-r)/(b+r); 
      mu = 3bmrdbpr*onembmr2inv
      p = bmrdbpr^2*onembpr2*onembmr2inv
      Piofk,Eofk,Em1mKdm = cel_bulirsch(k2inv,kc,p,1+mu,one(T),one(T),p+mu,k2c,zero(T))
      sqonembmr2 = sqrt(onembmr2)
#      Lambda1 = 2*sqonembmr2*(onembpr2*Piofk -(4-7r^2-b^2)*Eofk)/3
      Lambda1 = 2*sqonembmr2*(onembpr2*Piofk -(4-7r^2-b^2)*Eofk)*third
      s2_grad[1] = -4*r*sqonembmr2*Eofk
#      s2_grad[2] = -4*r/3*sqonembmr2*(Eofk - 2*Em1mKdm)
      s2_grad[2] = -4*r*third*sqonembmr2*(Eofk - 2*Em1mKdm)
    else
      # b+r = 1 or k^2=1, Case 4 (extending r up to 1)
#      Lambda1 = 2*acos(1.0-2.0*r)-4/3*(3+2r-8r^2)*sqrt(r*(1-r))-2pi*convert(T,r>0.5) 
      Lambda1 = 2*acos(1.0-2.0*r)-4*third*(3+2r-8r^2)*sqrt(r*(1-r))-2pi*convert(T,r>0.5) 
      s2_grad[1] = -8*r*sqrt(r*(1-r))
#      s2_grad[2] = -s2_grad[1]/3
      s2_grad[2] = -s2_grad[1]*third
      Eofk = one(T)
      Em1mKdm = one(T)
    end
  end
end
#flux = ((1.0-convert(T,r>b))*2pi-Lambda1)/3
flux = ((1.0-convert(T,r>b))*2pi-Lambda1)*third
return flux,Eofk,Em1mKdm
end
