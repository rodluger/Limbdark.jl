# This is the Kuncir (1961) adaptive simpson algorithm
# Modified to work with a vector function.
function simpson_vec(a::T, b::T, f::Function, I_of_f::Array{T,1}, i::T, epsilon::T, N::Int64, nf::Int64) where {T <: Real}

# value a, b, epsilon, N; 
# integer N;
# real a, b, I_of_f, i, epsilon; 
# real procedure f;
# comment This procedure integrates the function f(x) using a modified 
# Simpson's Rule quadrature formula. The quadrature is performed over j 
# subintervals of [a,b] forming the total area I_of_f. Convergence in each 
# subinterval of length (b-a)/2^n is indicated when the relative difference 
# between successive three-point and five-point area approximations
#  A_{3,j} = (b-a)(g_o + 4g_2 + g_4)/(3*2^{n+1})
#  A_{5,j} = (b-a)(g_o + 4g_1 + 2g_2 + 4g_3 + g_4)/(3*2{n+2})
# is less than or equal to an appropriate portion of the over-all tolerance 
# epsilon (i.e., |(A_{5,j} - A_{3,j})/A_{5,j}| \le  \epsilon/2^n with n \le N). 
# SIMPSON will reduce the size of each interval until this condition is satisfied.
# Complete integration over [a,b] is indicated by i = b. A value 
# a =< i < b indicates that the integration was terminated, leaving I_of_f the true 
# area under f in [a,i]. Further integration over [i,b] will necessitate either 
# the assignment of a larger N, a larger epsilon, or an integral substitution reducing 
# the slope of the integrand in that interval. It is recommended that this 
# procedure be used between known integrand maxima and minima.

third = one(T)/3
if VERSION >= v"0.7"
  g=Array{Float64,2}(undef,nf,5)
  A=Array{Float64,2}(undef,nf,3)
  S=Array{Float64,3}(undef,nf,3,N)
else
  g=Array{Float64}(nf,5)
  A=Array{Float64}(nf,3)
  S=Array{Float64}(nf,3,N)
end
f_of_x = zeros(T,nf)
I_of_f .= zero(T)
i=zero(T)
m=0; n=0
f(a,f_of_x)
@inbounds for j=1:nf
  g[j,1] = f_of_x[j]
end
f(0.5*(a+b),f_of_x)
@inbounds for j=1:nf
  g[j,3] = f_of_x[j]
end
f(b,f_of_x)
@inbounds for j=1:nf
  g[j,5] = f_of_x[j]
end
bma = b-a
@inbounds for j=1:nf
  A[j,1] = 0.5*bma*(g[j,1]+4*g[j,3]+g[j,5])
end
@label AA
d = 2^n
h = 0.25*bma/d

f(a+h*(4*m+1),f_of_x)
@inbounds for j=1:nf
  g[j,2] = f_of_x[j]
end
f(a+h*(4*m+3),f_of_x)
@inbounds for j=1:nf
  g[j,4] = f_of_x[j]
end
#println("g: ",g[1,:])
@inbounds for j=1:nf
  A[j,2] = h*(g[j,1]+4*g[j,2]+g[j,3])
  A[j,3] = h*(g[j,3]+4*g[j,4]+g[j,5])
end
maxdiff = zero(T)
@inbounds for j=1:nf
  maxa231 = abs((A[j,2]+A[j,3])-A[j,1])
  if maxa231 > maxdiff
    maxdiff = maxa231
  end
end
if maxdiff > 3*epsilon
  m *=2
  n +=1
  if n > N
    @goto CC
  end
  @inbounds for j=1:nf
    A[j,1] = A[j,2]
    S[j,1,n] = A[j,3]
    S[j,2,n] = g[j,4]
    S[j,3,n] = g[j,5]
    g[j,5] = g[j,3]
    g[j,3] = g[j,2]
  end
  @goto AA
else
  @inbounds for j=1:nf
    I_of_f[j] += (A[j,2]+A[j,3])*third
  end
#  println("I_of_f: ",I_of_f[1]," b-a: ",b-a," A: ",A," h: ",h," n: ",n)
  m += 1
  i = a+m*bma/d
  @label BB
  if m == 2*div(m,2)
    m = div(m,2)
    n -= 1
    @goto BB
  end
  if (m != 1) || (n != 0)
    @inbounds for j=1:nf
      A[j,1] = S[j,1,n]
      g[j,1] = g[j,5]
      g[j,3] = S[j,2,n]
      g[j,5] = S[j,3,n]
    end
    @goto AA
  end
end
@label CC
#println("I_of_f: ",I_of_f[1]," b-a: ",b-a)
return I_of_f
end
