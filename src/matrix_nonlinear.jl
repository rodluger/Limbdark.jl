# Computes matrix for finding coefficients of nonlinear
# limb-darkening in terms of polynomial limb-darkening:
function coeff_power2(alpha::Float64,N::Int64)
mat = zeros(N+1,N+1)
# n >= 2 terms
for i=3:N+1,k=3:N+1
  mat[i,k] = 4*(3-(i-1+k-1)+2*(i-1)*(k-1))/(i-1+k-1-3)/((i-1+k-1)^2-1)
end
# Cross terms n=0,1 with n >= 2
for i=1:N+1
  mat[i,1] = 2/(1-(i-1)^2)
  mat[1,i] = 2/(1-(i-1)^2)
end
# Now for n=0,1 terms:
mat[1,1] = 1.0
mat[1,2] = 0.5
mat[2,1] = 0.5
mat[2,2] = 1./3.

# Now, solve for power-2 limb-darkening coefficient:

# Take dot product of (1-mu^\alpha) to be multiplied by c:
imu_dot_g_n = zeros(N+1)
imu_dot_g_n[1] = alpha/(1+alpha)
imu_dot_g_n[2] = 0.5-1/(2+alpha)
for i=3:N+1
  n=i-1
  imu_dot_g_n[i] = 2*alpha*(1-alpha-n*(n+2))/(n^2-1)/((alpha+n)^2-1)
end
return \(mat,imu_dot_g_n)
end

# Now, loop over various alphas:
alpha = [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99,1.01,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8]
nalpha = length(alpha)
nord = 15
coeff = zeros(nalpha,nord,nord+1)
scat = zeros(nalpha,nord)
nmu = 1000; mu = linspace(1e-3,1.0,nmu)
using PyPlot
c = 0.5
for ialpha=1:nalpha
  fnon = 1-mu.^alpha[ialpha]
  for i=2:nord
    fn = zeros(i+1,nmu)
#    fn[1,:] .= big(1.0)
    fn[1,:] .= 1.0
    fn[2,:] .= mu
    for j=2:i
      fn[j+1,:] .= (j+2)*mu.^j-j*mu.^(j-2)
    end
    ci = coeff_power2(alpha[ialpha],i)
    fmod = zeros(Float64,nmu)
    for j=1:i+1
      fmod += ci[j]*fn[j,:]
    end
    clf()
    plot(mu,1-c*fnon)
    plot(mu,1-c*fmod)
    plot(mu,c*(fmod-fnon))
    scat[ialpha,i]=convert(Float64,std(fmod-fnon))
    coeff[ialpha,i,1:i+1]=convert(Array{Float64,1},ci)
    println("alpha: ",convert(Float64,alpha[ialpha])," i: ",i," scatter: ",scat[ialpha,i]," max: ",maximum(abs,convert(Array{Float64,1},fmod-fnon)))
  end
end
