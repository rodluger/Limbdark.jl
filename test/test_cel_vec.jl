# Tests the vectorized version of cel:
#include("cel_bulirsch.jl")
include("../src/cel_bulirsch.jl")
if VERSION >= v"0.7"
  using Test
else 
  using Base.Test
end

function test_cel!(kc,p,a,b)
k2 = 1.0-kc^2
ell1 = zeros(3)
nphi = 100000
if p < 0
  phi0 = asin(sqrt(1.0/(1.0-p)))
  phi1 = reverse(phi0 .-logarithmspace(log10(1e-8),log10(phi0),nphi+1))
  dphi = phi1[2:nphi]-phi1[1:nphi-1]
  phi = .5*(phi1[2:nphi]+phi1[1:nphi-1])
  phi1 = phi0 .+logarithmspace(log10(1e-8),log10(pi/2-phi0),nphi+1)
  dphi = [dphi;phi1[2:nphi]-phi1[1:nphi-1]]
  phi = [phi;.5*(phi1[2:nphi]+phi1[1:nphi-1])]
# Finally, evaluate numerically:
  ell3 = zeros(3); cphi2 = cos.(phi).^2; sphi2=sin.(phi).^2
  den = dphi./sqrt.(cphi2+kc*kc*sphi2)
  ell3[1] = sum((a[1]*cphi2+b[1]*sphi2)./(cphi2+p*sphi2).*den)
  ell3[2] = sum((a[2]*cphi2+b[2]*sphi2).*den)
  ell3[3] = sum((a[3]*cphi2+b[3]*sphi2).*den)
  p0 = sqrt((kc*kc-p)/(1-p))
  a0 = (a[1]-b[1])/(1-p)
  b0 = -(b[1]-a[1]*p)/(1-p)^2*(1-kc^2)/p0+a0*p0
  p = p0^2
  a[1] = a0
  b[1] = b0*p0
else
  dphi = .5*pi/nphi
  phi = linearspace(0.5*dphi,pi/2-.5*dphi,nphi)
  ell3 = zeros(3); cphi2 = cos.(phi).^2; sphi2=sin.(phi).^2
  den = dphi./sqrt.(cphi2+kc*kc*sphi2)
  ell3[1] = sum((a[1]*cphi2+b[1]*sphi2)./(cphi2+p*sphi2).*den)
  ell3[2] = sum((a[2]*cphi2+b[2]*sphi2).*den)
  ell3[3] = sum((a[3]*cphi2+b[3]*sphi2).*den)
end

# First, evaluate one by one:
ell1[1] = cel_bulirsch(k2,kc,p,a[1],b[1])
ell1[2] = cel_bulirsch(k2,kc,1.0,a[2],b[2])
ell1[3] = cel_bulirsch(k2,kc,1.0,a[3],b[3])
# Next, evaluate all three at once:
ell2 = cel_bulirsch(k2,kc,p,a[1],a[2],a[3],b[1],b[2],b[3])

#println("cel_bul: ",ell1)
#println("cel_alt: ",ell1_alt)
#println("numeric: ",ell3)
#println("sca-num: ",ell1-ell3)
return ell1,ell2,ell3
end

@testset "cel_bulirsch" begin
  # Test negative p values:
  p = -rand(); a=rand(3); b=rand(3)
  kc = rand()
  ell1,ell2,ell3 = test_cel!(kc,p,a,b)
  println("ell1: ",ell1)
  println("ell2: ",ell2)
  println("ell3: ",ell3)
  for i=1:3
    @test isapprox(ell1[i],ell2[i])
    @test isapprox(ell1[i],ell3[i],atol=1e-3)
  end
  # Test positive p values:
  p = rand(); a=rand(3); b=rand(3)
  kc = rand()
  ell1,ell2,ell3 = test_cel!(kc,p,a,b)
  println("ell1: ",ell1)
  println("ell2: ",ell2)
  println("ell3: ",ell3)
  for i=1:3
    @test isapprox(ell1[i],ell2[i])
    @test isapprox(ell1[i],ell3[i],atol=1e-3)
  end
  # Test small kc values:
  p = rand(); a=rand(3); b=rand(3)
  kc = 1e-4
  ell1,ell2,ell3 = test_cel!(kc,p,a,b)
  println("ell1: ",ell1)
  println("ell2: ",ell2)
  println("ell3: ",ell3)
  for i=1:3
    @test isapprox(ell1[i],ell2[i])
    @test isapprox(ell1[i],ell3[i],atol=1e-3)
  end
  # Test large kc values:
  p = rand(); a=rand(3); b=rand(3)
  kc = 1.0-1e-4
  ell1,ell2,ell3 = test_cel!(kc,p,a,b)
  println("ell1: ",ell1)
  println("ell2: ",ell2)
  println("ell3: ",ell3)
  for i=1:3
    @test isapprox(ell1[i],ell2[i])
    @test isapprox(ell1[i],ell3[i],atol=1e-3)
  end
end
