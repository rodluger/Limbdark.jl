using PyPlot
nu = 3
mu = linspace(0,1,10000)
ntrial = 100

function iofmu(u,mu)
  intensity = ones(mu)
  for i=1:length(u)
    intensity .-= u[i]*(1 .-mu).^i
  end
  return intensity
end

function diofmudmu(u,mu)
  didmu = zeros(mu)
  for i=1:length(u)
    didmu .+= i*u[i]*(1 .-mu).^(i-1)
  end
  return didmu
end

coeff = zeros(nu,ntrial)
u = zeros(nu)
clf()
for i=1:ntrial
  bad = true
  while bad
    u = -1.0+2.0*rand(nu)
    test = minimum(iofmu(u,mu)) < 0.0 || minimum(diofmudmu(u,mu)) <= 0.0
    if test
      bad = true
#      plot(mu,iofmu(u,mu),linestyle=":")
    else
      bad = false
      plot(mu,iofmu(u,mu))
      read(STDIN,Char)
      clf()
    end
  end
end
