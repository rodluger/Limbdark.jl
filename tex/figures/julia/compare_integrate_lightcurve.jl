# This is performance and accuracy tests for the integrate_lightcurve function

using BenchmarkTools
using PyPlot
using PyCall

# Python modules used for comparison
exo = pyimport("exoplanet")
tt = pyimport("theano.tensor")
starry = pyimport("starry")

include("../../../src/integrate_lightcurve.jl")
include("../../../test/loglinspace.jl")

println("\nRunning comparison...\n")

function run_ALFM_test(r::Ty, b0::Ty, aonr::Ty, period::Ty, u_1::Ty, u_2::Ty, trials::Int64) where {Ty <: Real}

  # TODO probably want to continue moving these out for sharing when running other code

  v = 2*pi*aonr/period # velocity
  trans_dur = 2*sqrt((1+r)^2-b0^2)/v*24 # approx transit duration in hours

  # Simulate transit data with a 2-minute cadence:
  tobs = 3*trans_dur  # Observation duration in hours
  nobs = convert(Int64,round(tobs * 30))  # Number of data points
  time = linearspace(-1.5*trans_dur ,1.5*trans_dur, nobs)  # Time array in units of hours


  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

    # Setup all the parameters
    trans = transit_init(r, b0, [u_1, u_2], false)  # Initialize transit structure
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
  end

  return nobs, timeTaken/trials
end

function run_exoplanet_test(orbit::PyObject, trials::Int64)

  u = tt.vector()
  b = tt.vector()
  r = tt.vector()
  #lc = StarryLightCurve(orbit=)


end


# Transit parameters. Some vary, some are fixed

# Ratio of semi-major axis to stellar radius:
# TODO may need to compute aonr?
aonr = 30.0
# Orbital period [d]:
period_vals = collect(logarithmspace(0.4, 4.0, 25)) # from about 2.5 to about 10000 days
r = 0.1
b = 0.3

# Quadratic limb darkening parameters
u_1 = 0.2
u_2 = 0.2


# Comparison configuration
runs = length(period_vals)
trialsPerConfig = 10 # number of trials to run for a given parameter configuration

# Structures for holding results
num_datapoints = zeros(runs)
limbdark_time = zeros(runs)
exoplanet_time = zeros(runs)

# Run the timing test
for i = 1:runs
  period = period_vals[mod(i-1, length(period_vals)) + 1]

  # with these parameters, it will assume solar mass and solar radius
  # TODO that probably isn't consistent with what we did above
  orbit = exo.orbits.KeplerianOrbit(period=period, b=b)

  # the number of data points is computed based on the parameters provided
  # this is the most interesting thing to graph
  # TODO use orbit data. mine is probably not a real orbit
  num_datapoints[i], limbdark_time[i] = run_ALFM_test(r, b, aonr, period, u_1, u_2, trialsPerConfig)

  #exoplanet_time[i] = run_exoplanet_test(orbit, trialsPerConfig)

end

# Plot Results
fig, axes = subplots(3,1, figsize=(8,12))

ax = axes[1]
ax.plot(num_datapoints, limbdark_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")
ax.set_title("r=0.1; b=0.3")
ax.set_xlabel("Data Points")
ax.set_ylabel("Time [s]")
ax.set_xscale(:log) # scales lineary to number of data points
ax.set_yscale(:log) # scales lineary to number of data points
ax.legend(loc="upper left",fontsize=10)

#=
ax = axes[2]
ax.plot(b_vals, limbdark_b_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")
ax.set_title("r=0.1; period=50d")
ax.set_xlabel("Impact Parameter")
ax.set_ylabel("Time [s]")
ax.legend(loc="upper left",fontsize=10)
=#
