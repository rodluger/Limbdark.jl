# This is performance and accuracy tests for the integrate_lightcurve function

using BenchmarkTools

if VERSION >= v"0.7"
  using Statistics
end

include("../../src/integrate_lightcurve.jl")
include("../loglinspace.jl")

@testset "Compare Performance and Accuracy of Integrate Lightcurve" begin

  function run_timing_test()

    # Specify parameters appropriate for KOI-984.01.
    # Plucking these from NexSci Kepler candidate archive.

    # Radius ratio:
    r = 0.025717
    # Ratio of semi-major axis to stellar radius:
    aonr = 19.94
    # Orbital period [d]:
    period = 4.28746668
    # Impact parameter:
    b0 = 0.374
    # Transit duration [hr] (this is time between first and fourth contact):
    #T0 = 1.5693 (this is value from catalog)

    # Check that these agree.  Compute velocity:
    v = 2*pi*aonr/period
    # Compute transit duration [convert to hr], approximate:
    T = 2*sqrt((1+r)^2-b0^2)/v*24
    #@test T ≈ T0 atol=0.001

    # Okay, now to simulate transit.  Assume a 2-minute cadence:
    tobs = 3*T  # Observation duration in hours
    nobs = convert(Int64,round(tobs * 30))  # Number of data points
    time = linearspace(-1.5*T ,1.5*T, nobs)  # Time array in units of hours

    u_n = [0.1914,0.5167]  # Quadratic limb-darkening parameters

    # Setup all the parameters
    trans = transit_init(r, b0, u_n, false)  # Initialize transit structure
    param = [0.0, v/24, b0]   # parameters of the transit: [t_0,v,b_0]
    t = Array{Float64,1}(undef,nobs)
    t .= time
    dt = 2/60  # Transit exposure time is 2 minutes, convert to hours.
    favg1 = Array{Float64,2}(undef, 7, nobs)
    nt = nobs
    tol = 1e-6
    maxdepth = 6
    neval_t = Array{Int64,1}(undef, nobs)
    depthmax = Array{Int64,1}(undef, nobs)

    integrate_lightcurve!(trans, param, t, dt, favg1, nt, tol, maxdepth, neval_t, depthmax)

    #println("Depthmax: $depthmax")
    #println("neval_t: $neval_t")
    #println("Depthmax: $depthmax")

  end

  # Run the timing test many times and benchmark the total amount of time
  expression, timeTaken, bytes, gctime, memallocs = @timed for i = 1:1000
  #@time for i = 1:1000
    run_timing_test()
  end

  println("Time taken: $timeTaken")
end
