# This is performance and accuracy tests for the integrate_lightcurve function

using BenchmarkTools
using PyPlot

if VERSION >= v"0.7"
  using Statistics
  using DelimitedFiles
end

include("../../../src/integrate_lightcurve.jl")
include("../../../test/loglinspace.jl")

println("\nRunning comparison...\n")

# TODO read in data from starry or exoplanet for comparison

function run_timing_test(r::Ty, b0::Ty, aonr::Ty, period::Ty, u_1::Ty, u_2::Ty, trials::Int64) where {Ty <: Real}

  # TODO probably want to continue moving these out for sharing when running other code

  v = 2*pi*aonr/period # velocity
  trans_dur = 2*sqrt((1+r)^2-b0^2)/v*24 # approx transit duration in hours

  # Simulate transit data with a 2-minute cadence:
  tobs = 3*trans_dur  # Observation duration in hours
  nobs = convert(Int64,round(tobs * 30))  # Number of data points
  time = linearspace(-1.5*trans_dur ,1.5*trans_dur, nobs)  # Time array in units of hours

  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

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
  end

  return nobs, timeTaken/trials
end


# Transit parameters. Some vary, some are fixed

# Ratio of semi-major axis to stellar radius:
aonr = 30.0
# Orbital period [d]:
period_vals = collect(linearspace(1, 1000.0, 20)) # = 4.28746668
#r = 0.025717 # radius ratio
r = 0.1
b_vals = [0.3]
#b_vals = collect(linearspace(0.1, 0.4, 4)) # impact parameter
#b[i+1] = sqrt(((i-500.0)/500.0*(1.0+2.0*r))^2) # use i-1
# Quadratic limb darkening parameters
u_1 = 0.2
u_2 = 0.2

# Comparison configuration
runs = length(period_vals) * length(b_vals) # number of parameter configurations we're trying
trialsPerConfig = 10 # number of trials to run for a given parameter configuration

# Structures for holding results
limbdark_time = zeros(runs)
num_datapoints = zeros(runs)
#results_aveflux = zeros(runs)

# Run the timing test
for i = 1:runs
  b = b_vals[mod(i-1, length(b_vals)) + 1]
  period = period_vals[mod(i-1, length(period_vals)) + 1]

  # the number of data points is computed based on the parameters provided
  # this is the most interesting thing to graph
  num_datapoints[i], limbdark_time[i] = run_timing_test(r, b, aonr,  period, trialsPerConfig)
end



# Plot Results
fig, axes = subplots(3,1, figsize=(8,12))

ax = axes[1]
ax.plot(num_datapoints, limbdark_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")

ax.set_title("r=0.1; b=0.3")
ax.set_xlabel("Data Points")
ax.set_xlabel("Data Points")

#ax.set_xscale(:log) # scales lineary to number of data points
ax.legend(loc="upper left",fontsize=10)

#=
ax = axes[2]
ax.plot(period_vals, limbdark_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")

ax.set_xlabel("Data Points")
ax.legend(loc="upper right",fontsize=10)
=#

# TODO plot other's data for comparison
