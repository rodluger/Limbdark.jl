# Tests the vectorized version of cel:
#include("cel_bulirsch.jl")
# include("../src/cel_bulirsch.jl")
#if VERSION >= v"0.7"
#  using Test
#else
#  using Base.Test
#end

import Limbdark: cel_bulirsch

# Randomizer seed
using Random
Random.seed!(42)

# Function which tests cel_bulirsch both numerically, one-by-one and vectorized:
function test_cel!(kc,p,a,b)
k2 = 1.0-kc^2
# Evaluate all three at once:
ell2 = cel_bulirsch(k2,kc,p,a[1],a[2],a[3],b[1],b[2],b[3])
ell1 = zeros(3)
nphi = 100000
if p < 0
#  Make same transformations as in cel_bulirsch for negative p:
  p0 = sqrt((kc*kc-p)/(1-p))
  a0 = (a[1]-b[1])/(1-p)
  b0 = -(b[1]-a[1]*p)/(1-p)^2*(1-kc^2)/p0+a0*p0
  a[1] = a0
  b[1] = b0*p0
  p = p0*p0
end
# Numerically compute the integrals:
dphi = .5*pi/nphi
phi = linearspace(0.5*dphi,pi/2-.5*dphi,nphi)
ell3 = zeros(3); cphi2 = cos.(phi).^2; sphi2=sin.(phi).^2
den = dphi./sqrt.(cphi2+kc*kc*sphi2)
ell3[1] = sum((a[1]*cphi2+b[1]*sphi2)./(cphi2+p*sphi2).*den)
ell3[2] = sum((a[2]*cphi2+b[2]*sphi2).*den)
ell3[3] = sum((a[3]*cphi2+b[3]*sphi2).*den)

# First, evaluate one by one:
ell1[1] = cel_bulirsch(k2,kc,p,a[1],b[1])
ell1[2] = cel_bulirsch(k2,kc,1.0,a[2],b[2])
ell1[3] = cel_bulirsch(k2,kc,1.0,a[3],b[3])


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
#    @test isapprox(ell1[i],ell3[i],atol=1e-6)
    @test isapprox(ell1[i],ell3[i])
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
#    @test isapprox(ell1[i],ell3[i],atol=1e-6)
    @test isapprox(ell1[i],ell3[i])
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
#    @test isapprox(ell1[i],ell3[i],atol=1e-6)
    @test isapprox(ell1[i],ell3[i])
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
#    @test isapprox(ell1[i],ell3[i],atol=1e-6)
    @test isapprox(ell1[i],ell3[i])
  end
end
