# Can we approximate non-linear limb-darkening with
# polynomial?

include("/Users/ericagol/Computer/Julia/regress.jl")

nmu = 1000; mu = linspace(0,1,nmu)
fnon = sqrt.(mu).^3-mu
using PyPlot
for i=1:20
  fn = zeros(i,nmu)
  fn[1,:] = 1.0
  for j=2:i
    fn[j,:] =mu.^(j-1)
  end
  coeff,cov = regress(fn,fnon,ones(nmu))
  fmod = zeros(nmu)
  for j=1:i
    fmod += coeff[j]*fn[j,:]
  end
  clf()
  plot(mu,fnon+mu)
  plot(mu,fmod+mu)
  plot(mu,fmod-fnon)
  println("i: ",i," scatter: ",std(fmod-fnon)," max: ",maximum(abs,fmod-fnon))
  println("minimum: ",minimum(fmod+mu)," maximum: ",maximum(fmod+mu))
  read(STDIN,Char)
end
